<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example RCAL-1234: <Fix a bug> -->
Resolves [RCAL-nnnn](https://jira.stsci.edu/browse/RCAL-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

<!-- if you can't perform these due to permissions, please ask a maintainer to do them -->
## Tasks
- [ ] **request a review from someone specific**, to avoid making the maintainers review every PR
- [ ] add a build milestone, i.e. `24Q4_B15` (use the [latest build](https://github.com/spacetelescope/romancal/milestones) if not sure)
- [ ] Does this PR change user-facing code / API? (if not, label with `no-changelog-entry-needed`)
  - [ ] write news fragment(s) in `changes/`: `echo "changed something" > changes/<PR#>.<changetype>.rst` (see below for change types)
  - [ ] update or add relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
  - [ ] [start a regression test](https://github.com/spacetelescope/RegressionTests/actions/workflows/romancal.yml) and include a link to the running job ([click here for instructions](https://github.com/spacetelescope/RegressionTests/blob/main/docs/running_regression_tests.md))
    - [ ] Do truth files need to be updated ("okified")?
      - [ ] **after the reviewer has approved these changes**, run `okify_regtests` to update the truth files
- [ ] if a JIRA ticket exists, [make sure it is resolved properly](https://github.com/spacetelescope/romancal/wiki/How-to-resolve-JIRA-issues)

<details><summary>news fragment change types...</summary>

  - ``changes/<PR#>.general.rst``: infrastructure or miscellaneous change
  - ``changes/<PR#>.docs.rst``
  - ``changes/<PR#>.stpipe.rst``
  - ``changes/<PR#>.associations.rst``
  - ``changes/<PR#>.scripts.rst``
  - ``changes/<PR#>.mosaic_pipeline.rst``
  - ``changes/<PR#>.skycell.rst``

  ## steps
  - ``changes/<PR#>.dq_init.rst``
  - ``changes/<PR#>.saturation.rst``
  - ``changes/<PR#>.refpix.rst``
  - ``changes/<PR#>.wfi18_transient.rst``
  - ``changes/<PR#>.linearity.rst``
  - ``changes/<PR#>.dark_current.rst``
  - ``changes/<PR#>.jump_detection.rst``
  - ``changes/<PR#>.ramp_fitting.rst``
  - ``changes/<PR#>.assign_wcs.rst``
  - ``changes/<PR#>.flatfield.rst``
  - ``changes/<PR#>.photom.rst``
  - ``changes/<PR#>.flux.rst``
  - ``changes/<PR#>.source_detection.rst``
  - ``changes/<PR#>.tweakreg.rst``
  - ``changes/<PR#>.skymatch.rst``
  - ``changes/<PR#>.outlier_detection.rst``
  - ``changes/<PR#>.resample.rst``
  - ``changes/<PR#>.source_catalog.rst``
</details>
