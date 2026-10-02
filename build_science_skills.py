#!/usr/bin/env python3
"""Build / refresh the Google DeepMind Science Skills vendored in .claude/skills/.

Upstream: https://github.com/google-deepmind/science-skills (Apache-2.0 code,
CC-BY-4.0 docs). The collection is vendored at a pinned release so every Claude
Code session on this repo gets the same 40 skills with no network step.

Usage (from the repo root, Python >= 3.11, stdlib only):

    python build_science_skills.py              # index + smoke-test what is vendored
    python build_science_skills.py --update     # re-sync from upstream at PINNED tag
    python build_science_skills.py --update --tag v1.3.0   # bump the pin
    python build_science_skills.py --no-smoke   # index only (offline)

Outputs:
    .claude/skills/<skill>/...                 the skills themselves
    .claude/skills/SCIENCE_SKILLS_INDEX.md     one-line index (name -> description)
    docs/SCIENCE_SKILLS_BUILD_REPORT.md        per-script smoke-test result

The smoke test runs every helper script with `--help` under `uv run`, which
resolves and caches each script's PEP 723 / pyproject dependencies. A script
that aborts only because an API key is absent is reported as CRED, not FAIL.
ASCII-only output (Windows cp1252 console).
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

UPSTREAM_REPO = "https://github.com/google-deepmind/science-skills.git"
PINNED_TAG = "v1.2.1"
PINNED_COMMIT = "68832757cbbf941c620b71df5756cf6e5cc287b0"

ROOT = Path(__file__).resolve().parent
SKILLS_DIR = ROOT / ".claude" / "skills"
INDEX_MD = SKILLS_DIR / "SCIENCE_SKILLS_INDEX.md"
REPORT_MD = ROOT / "docs" / "SCIENCE_SKILLS_BUILD_REPORT.md"
LICENSE_DST = SKILLS_DIR / "LICENSE-google-deepmind-science-skills.txt"
SKILL_LICENSES_DST = SKILLS_DIR / "SKILL_LICENSES.md"

USAGE_PATTERN = re.compile(r"^\s*usage:", re.I | re.M)
CRED_PATTERN = re.compile(r"(API_KEY|USER_EMAIL)[^\n]*not set|Missing credential", re.I)
JUNK_DIRS = {".venv", "__pycache__", ".pytest_cache", ".mypy_cache"}
JUNK_FILES = {"uv.lock"}


def log(msg: str) -> None:
    print(msg.encode("ascii", "replace").decode("ascii"), flush=True)


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, text=True, capture_output=True, **kw)


# --------------------------------------------------------------------------- sync

def scrub(path: Path) -> None:
    """Remove build artefacts so they never land in git."""
    for p in sorted(path.rglob("*"), key=lambda q: -len(q.parts)):
        if p.is_dir() and p.name in JUNK_DIRS:
            shutil.rmtree(p, ignore_errors=True)
        elif p.is_file() and (p.name in JUNK_FILES or p.suffix == ".pyc"):
            p.unlink(missing_ok=True)


def update_from_upstream(tag: str, expected_commit: str | None) -> str:
    if shutil.which("git") is None:
        sys.exit("git is required for --update")
    with tempfile.TemporaryDirectory(prefix="science-skills-") as tmp:
        clone = Path(tmp) / "upstream"
        log(f"[update] cloning {UPSTREAM_REPO} @ {tag}")
        r = run(["git", "clone", "--quiet", "--depth", "1", "--branch", tag, UPSTREAM_REPO, str(clone)])
        if r.returncode:
            sys.exit(f"git clone failed:\n{r.stderr}")
        commit = run(["git", "-C", str(clone), "rev-parse", "HEAD"]).stdout.strip()
        if expected_commit and commit != expected_commit:
            sys.exit(f"[update] commit mismatch for {tag}: got {commit}, pinned {expected_commit}. "
                     "Pass --commit to accept the new hash.")
        src_skills = clone / "skills"
        upstream_names = sorted(p.name for p in src_skills.iterdir() if p.is_dir())
        SKILLS_DIR.mkdir(parents=True, exist_ok=True)
        for name in upstream_names:
            dst = SKILLS_DIR / name
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src_skills / name, dst)
        shutil.copy2(clone / "LICENSE", LICENSE_DST)
        shutil.copy2(clone / "SKILL_LICENSES.md", SKILL_LICENSES_DST)
        scrub(SKILLS_DIR)
        log(f"[update] synced {len(upstream_names)} skills at {tag} ({commit[:12]})")
        return commit


# --------------------------------------------------------------------------- index

def read_frontmatter(skill_md: Path) -> tuple[str, str]:
    text = skill_md.read_text(encoding="utf-8")
    m = re.match(r"^---\n(.*?)\n---", text, re.S)
    if not m:
        return skill_md.parent.name, ""
    fm = m.group(1)
    name = re.search(r"^name:\s*(.+)$", fm, re.M)
    desc_lines: list[str] = []
    capture = False
    for line in fm.splitlines():
        if re.match(r"^description:\s*(>-?|\|)?\s*(.*)$", line):
            capture = True
            tail = re.match(r"^description:\s*(?:>-?|\|)?\s*(.*)$", line).group(1).strip()
            if tail:
                desc_lines.append(tail)
            continue
        if capture:
            if re.match(r"^[A-Za-z_-]+:", line):
                break
            desc_lines.append(line.strip())
    desc = " ".join(x for x in desc_lines if x).strip()
    return (name.group(1).strip() if name else skill_md.parent.name), desc


def vendored_skills() -> list[Path]:
    return sorted(p for p in SKILLS_DIR.iterdir() if p.is_dir() and (p / "SKILL.md").exists())


def write_index(tag: str, commit: str) -> int:
    rows = []
    for d in vendored_skills():
        name, desc = read_frontmatter(d / "SKILL.md")
        scripts = sorted((d / "scripts").glob("*.py")) if (d / "scripts").exists() else []
        rows.append((name, d.name, len(scripts), desc))
    lines = [
        "# Google DeepMind Science Skills -- vendored index",
        "",
        f"Upstream: {UPSTREAM_REPO}  ",
        f"Pinned: `{tag}` (`{commit[:12]}`)  ",
        f"Generated: {dt.date.today().isoformat()} by `build_science_skills.py` -- do not edit by hand.",
        "",
        "| # | Skill (frontmatter name) | Directory | Scripts | Description |",
        "|---|---|---|---|---|",
    ]
    for i, (name, dname, n, desc) in enumerate(rows, 1):
        desc = desc.replace("|", "\\|")
        lines.append(f"| {i} | `{name}` | `{dname}/` | {n} | {desc} |")
    INDEX_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"[index] {len(rows)} skills -> {INDEX_MD.relative_to(ROOT)}")
    return len(rows)


# --------------------------------------------------------------------------- smoke

def smoke_test(timeout: int) -> list[tuple[str, str, str, int, str]]:
    uv = shutil.which("uv") or str(Path.home() / ".local" / "bin" / "uv")
    if not Path(uv).exists():
        sys.exit("uv not found. Install: curl -LsSf https://astral.sh/uv/install.sh | sh "
                 "(see .claude/skills/uv/SKILL.md), or re-run with --no-smoke.")
    results = []
    for d in vendored_skills():
        scripts_dir = d / "scripts"
        if not scripts_dir.exists():
            continue
        project = (d / "pyproject.toml").exists()
        for script in sorted(scripts_dir.glob("*.py")):
            cmd = [uv, "run"] + (["--project", str(d)] if project else []) + [str(script), "--help"]
            t0 = dt.datetime.now()
            try:
                r = run(cmd, cwd=str(d), timeout=timeout)
                rc, out = r.returncode, (r.stderr or "") + (r.stdout or "")
            except subprocess.TimeoutExpired:
                rc, out = 124, "timeout"
            secs = int((dt.datetime.now() - t0).total_seconds())
            if rc == 0 or USAGE_PATTERN.search(out):
                status = "PASS"  # printed usage (pubmed_api.py exits 1 on --help by design)
            elif CRED_PATTERN.search(out):
                status = "CRED"
            else:
                status = "FAIL"
            last = next((ln for ln in reversed(out.strip().splitlines()) if ln.strip()), "")
            rel = str(script.relative_to(SKILLS_DIR)).replace(os.sep, "/")
            results.append((rel, "project" if project else "script", status, secs, last.strip()[:160]))
            log(f"[smoke] {status:4s} {secs:4d}s {rel}")
    scrub(SKILLS_DIR)
    return results


def write_report(results, tag: str, commit: str, n_skills: int) -> None:
    counts = {k: sum(1 for r in results if r[2] == k) for k in ("PASS", "CRED", "FAIL")}
    lines = [
        "# Science Skills build report",
        "",
        f"Generated {dt.datetime.now().strftime('%Y-%m-%d %H:%M')} by `build_science_skills.py`.  ",
        f"Upstream `{tag}` (`{commit[:12]}`); {n_skills} skills vendored in `.claude/skills/`.",
        "",
        "Each helper script was executed as `uv run [--project <skill>] <script> --help`.",
        "PASS = resolved dependencies and printed usage. CRED = aborted only because an API key",
        "or e-mail is not in `~/.env` (add it per `.claude/skills/credentials/SKILL.md`).",
        "FAIL = anything else; see the last line captured.",
        "",
        f"| PASS | CRED | FAIL | total |",
        f"|---|---|---|---|",
        f"| {counts['PASS']} | {counts['CRED']} | {counts['FAIL']} | {len(results)} |",
        "",
        "| Script | Mode | Status | s | Last line |",
        "|---|---|---|---|---|",
    ]
    for rel, mode, status, secs, last in results:
        last_safe = last.replace("|", "\\|")
        lines.append(f"| `{rel}` | {mode} | {status} | {secs} | {last_safe} |")
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"[report] PASS={counts['PASS']} CRED={counts['CRED']} FAIL={counts['FAIL']} -> {REPORT_MD.relative_to(ROOT)}")


# --------------------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="re-sync skills from upstream at --tag")
    ap.add_argument("--tag", default=PINNED_TAG, help=f"upstream release tag (default {PINNED_TAG})")
    ap.add_argument("--commit", default=None, help="expected commit for --tag (default: pinned hash when tag is the pinned tag)")
    ap.add_argument("--no-smoke", action="store_true", help="skip the uv smoke test")
    ap.add_argument("--timeout", type=int, default=900, help="per-script timeout in seconds")
    a = ap.parse_args()

    commit = PINNED_COMMIT
    if a.update:
        expected = a.commit or (PINNED_COMMIT if a.tag == PINNED_TAG else None)
        commit = update_from_upstream(a.tag, expected)
    if not SKILLS_DIR.exists():
        sys.exit(f"{SKILLS_DIR} missing -- run with --update first")
    n = write_index(a.tag, commit)
    if a.no_smoke:
        return 0
    results = smoke_test(a.timeout)
    write_report(results, a.tag, commit, n)
    return 1 if any(r[2] == "FAIL" for r in results) else 0


if __name__ == "__main__":
    sys.exit(main())
