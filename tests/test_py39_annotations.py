import ast
from pathlib import Path

FUTURE_IMPORT = "annotations"
SEARCH_ROOTS = ("src", "tests")


def _annotation_uses_pep604_union(annotation: ast.AST) -> bool:
    return any(
        isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr)
        for node in ast.walk(annotation)
    )


def _analyze_annotations(path: Path) -> tuple[bool, bool]:
    tree = ast.parse(path.read_text())
    has_future_import = any(
        isinstance(node, ast.ImportFrom)
        and node.module == "__future__"
        and any(alias.name == FUTURE_IMPORT for alias in node.names)
        for node in tree.body
    )
    for node in ast.walk(tree):
        if isinstance(node, ast.arg) and node.annotation:
            if _annotation_uses_pep604_union(node.annotation):
                return True, has_future_import
        elif isinstance(node, ast.AnnAssign):
            if _annotation_uses_pep604_union(node.annotation):
                return True, has_future_import
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.returns:
            if _annotation_uses_pep604_union(node.returns):
                return True, has_future_import
    return False, has_future_import


def test_pep604_annotations_are_safe_on_python39():
    root = Path(__file__).resolve().parents[1]
    checked = []
    missing = []
    for search_root in SEARCH_ROOTS:
        for path in sorted((root / search_root).rglob("*.py")):
            has_pep604_annotations, has_future_import = _analyze_annotations(path)
            if not has_pep604_annotations:
                continue
            checked.append(path.relative_to(root).as_posix())
            if not has_future_import:
                missing.append(path.relative_to(root).as_posix())

    assert checked
    assert not missing
