# Contributing

Contributions should preserve the project's narrow scientific scope and mathematical transparency.

A useful change should normally include:

- a clear statement of the mathematical or user-facing behavior
- tests for new public behavior
- documentation for public API changes
- validation when numerical accuracy or performance claims change
- no generated caches, build artifacts, benchmark bundles, or local backup directories in commits.

For larger numerical changes, compare against an analytic field whenever possible. Avoid weakening a mathematical test merely to make a new implementation pass. First determine whether the implementation or the test encodes the wrong mathematics.
