# CLAUDE.md

Project-level instructions for any future Claude Code session working in this
repository. Read this before editing code.

## Versioning workflow (MANDATORY)

This project tracks its version in two places that must stay in lock-step:

1. The `version` field in `Cargo.toml` (`[package]` block).
2. The top of `CHANGELOG.md`.

**Every time you make a change to the codebase**, before reporting the work as
done, you MUST:

1. **Increment the patch number** (the last component) of the version in
   `Cargo.toml` by exactly one.
   - Example: `0.9.1` → `0.9.2`, `0.9.2` → `0.9.3`, etc.
   - Do not bump minor or major numbers on your own — that is a human decision.
   - Do this for any source change: bug fix, refactor, feature, doc edit that
     ships with the code, build/config change. If the only thing you touched
     is `CHANGELOG.md` itself, do not bump again.

2. **Add a new entry to `CHANGELOG.md`** under a new version heading matching
   the new `Cargo.toml` version. Use the
   [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format:

   ```markdown
   ## [X.Y.Z] - YYYY-MM-DD

   ### Added
   - ...

   ### Changed
   - ...

   ### Fixed
   - ...

   ### Removed
   - ...
   ```

   Only include the subsections that apply. Use today's date in ISO format.
   Write entries from the user's perspective — what changed about the program,
   not which file you edited.

3. Commit `Cargo.toml`, `Cargo.lock` (if regenerated), and `CHANGELOG.md`
   together with the rest of the change. They are one logical unit.

If you forget to bump or update the changelog, the change is not complete.

## Files that are tracked by git despite the `*.md` ignore

`.gitignore` ignores all markdown files except an explicit allowlist. The
currently allowlisted files are:

- `README.md`
- `results_format.md`
- `CHANGELOG.md`
- `CLAUDE.md`

If you add a new markdown file that should be tracked, add it to the allowlist
in `.gitignore`. Otherwise it will be silently ignored.

Do not make any commits or push to github on your own. This is currently a human task.

## Open work

See `version1_todo.md` for the road-to-1.0 punch list. That file is
intentionally untracked (covered by the `*.md` ignore) and lives only locally.
