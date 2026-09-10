<!-- If this PR addresses a JIRA ticket: -->
<!-- Resolves [RCAL-nnnn](https://jira.stsci.edu/browse/RCAL-nnnn) -->

<!-- If this PR will close an existing GitHub issue (that is not already attached to a JIRA ticket): -->
<!-- Closes # -->

<!-- Describe your changes here: -->

## Description

This change ...

<!-- If you can't perform these tasks due to permissions, reach out to a maintainer. -->
## Tasks

- [ ] **request a review from someone specific**, to avoid making the maintainers review every PR
- [ ] add a build milestone, i.e. `24Q4_B15` (use the [latest build](https://github.com/spacetelescope/romancal/milestones) if not sure)
- [ ] Does this PR change user-facing code / API? (if not, label with `no-changelog-entry-needed`)
  - [ ] write news fragment(s) in `changes/`: `echo "changed something" > changes/<PR#>.<changetype>.rst` (see [changelog readme](https://github.com/spacetelescope/romancal/blob/main/changes/README.rst) for instructions)
    - if your change breaks existing functionality, also add a `changes/<PR#>.breaking.rst` news fragment
  - [ ] update or add relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
  - [ ] [start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/romancal.yml) and include a link to the running job ([click here for instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md))
    - [ ] Do truth files need to be updated ("okified")?
      - [ ] **after the reviewer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/romancal/wiki/How-to-resolve-JIRA-issues)

## Generative AI Usage Disclosure

<!-- If generative AI or LLMs were used in the process of making this change, describe their use here. -->
<!-- Otherwise, indicate "No genAI tools used". -->
