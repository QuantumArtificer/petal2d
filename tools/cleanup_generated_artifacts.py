#!/usr/bin/env python3
"""Remove generated PETAL2D development artifacts without deleting tracked files.

By default this script performs a dry run. Pass ``--apply`` to delete candidates.
Any candidate containing a Git-tracked file is skipped and reported for manual
review.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def tracked_paths(root: Path) -> set[Path]:
    try:
        out = subprocess.check_output(
            ["git", "ls-files", "-z"], cwd=root, stderr=subprocess.DEVNULL
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return set()
    return {root / p.decode() for p in out.split(b"\0") if p}


def contains_tracked(path: Path, tracked: set[Path]) -> bool:
    if path.is_file() or path.is_symlink():
        return path in tracked
    prefix = str(path.resolve()) + "/"
    return any(str(p.resolve()).startswith(prefix) for p in tracked)


def candidates(root: Path) -> list[Path]:
    out: set[Path] = set()
    explicit = [
        root / "build",
        root / "dist",
        root / ".pytest_cache",
        root / "docs" / "_build",
        root / "validation" / "results",
        root / "benchmarks" / "results",
        root / "src" / "gasp2d",
        root / "src" / "gasp2d.egg-info",
    ]
    out.update(p for p in explicit if p.exists())

    for pattern in [
        ".petal2d_backup_*",
        "petal2d_*analysis*",
        "petal2d_*analysis*.zip",
        "petal2d_*diagnostic*.json",
        "petal2d_*diagnostic*.py",
        "package_petal2d_*results.sh",
        "*.egg-info",
        "src/*.egg-info",
    ]:
        out.update(root.glob(pattern))

    out.update(root.rglob("__pycache__"))
    out.update(root.rglob("*.pyc"))
    out.update(root.rglob("*.pyo"))
    return sorted(out, key=lambda p: (len(p.parts), str(p)), reverse=True)


def remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="Repository root")
    parser.add_argument(
        "--apply", action="store_true", help="Actually delete untracked candidates"
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    tracked = tracked_paths(root)
    mode = "APPLY" if args.apply else "DRY RUN"
    print(f"PETAL2D cleanup: {mode}")
    print(f"Repository: {root}\n")

    removed = skipped = 0
    for path in candidates(root):
        # A parent may already have been removed.
        if not path.exists() and not path.is_symlink():
            continue
        rel = path.relative_to(root)
        if contains_tracked(path, tracked):
            print(f"SKIP tracked content: {rel}")
            skipped += 1
            continue
        if args.apply:
            remove(path)
            print(f"REMOVED: {rel}")
        else:
            print(f"WOULD REMOVE: {rel}")
        removed += 1

    print(f"\nCandidates {'removed' if args.apply else 'identified'}: {removed}")
    print(f"Skipped because Git-tracked content was present: {skipped}")
    if not args.apply:
        print("\nRun again with --apply after reviewing the list.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
