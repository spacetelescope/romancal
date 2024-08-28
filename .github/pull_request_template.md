<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example RCAL-1234: <Fix a bug> -->
Resolves [RCAL-nnnn](https://jira.stsci.edu/browse/RCAL-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

**Checklist**
- [ ] for a public change, added a towncrier news fragment in `changes/` <details><summary>`echo "changed something" > changes/<PR#>.<changetype>.rst`</summary>

    - ``changes/<PR#>.general.rst``: infrastructure or miscellaneous change
    - ``changes/<PR#>.docs.rst``
    - ``changes/<PR#>.stpipe.rst``
    - ``changes/<PR#>.associations.rst``
    - ``changes/<PR#>.scripts.rst``
    - ``changes/<PR#>.mosaic_pipeline.rst``
    - ``changes/<PR#>.patch_match.rst``

    ## steps
    - ``changes/<PR#>.dq_init.rst``
    - ``changes/<PR#>.saturation.rst``
    - ``changes/<PR#>.refpix.rst``
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
- [ ] updated relevant tests
- [ ] updated relevant documentation
- [ ] updated relevant milestone(s)
- [ ] added relevant label(s)
- [ ] ran regression tests, post a link to the Jenkins job below. [How to run regression tests on a PR](https://github.com/spacetelescope/romancal/wiki/Running-Regression-Tests-Against-PR-Branches)
