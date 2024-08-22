<!-- If this PR closes a JIRA ticket, make sure the title starts with the JIRA issue number,
for example RCAL-1234: <Fix a bug> -->
Resolves [RCAL-nnnn](https://jira.stsci.edu/browse/RCAL-nnnn)

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #

<!-- describe the changes comprising this PR here -->
This PR addresses ...

**Checklist**
- [ ] for a public change, added a [towncrier news fragment](https://towncrier.readthedocs.io/en/stable/tutorial.html#creating-news-fragments) <details><summary>`changes/<PR#>.<changetype>.rst`</summary>

    - ``changes/<PR#>.docs.rst``: documentation change
    - ``changes/<PR#>.general.rst``: infrastructure or miscellaneous change
    - ``changes/<PR#>.scripts.rst``: change to scripts
    - ``changes/<PR#>.stpipe.rst``: change to `stpipe` interface

    ## steps
    - ``changes/<PR#>.associations.rst``
    - ``changes/<PR#>.dark_current.rst``
    - ``changes/<PR#>.dq_init.rst``
    - ``changes/<PR#>.flux.rst``
    - ``changes/<PR#>.flatfield.rst``
    - ``changes/<PR#>.jump_detection.rst``
    - ``changes/<PR#>.linearity.rst``
    - ``changes/<PR#>.mosaic_pipeline.rst``
    - ``changes/<PR#>.outlier_detection.rst``
    - ``changes/<PR#>.patch_match.rst``
    - ``changes/<PR#>.photom.rst``
    - ``changes/<PR#>.ramp_fitting.rst``
    - ``changes/<PR#>.refpix.rst``
    - ``changes/<PR#>.refpix.rst``
    - ``changes/<PR#>.resample.rst``
    - ``changes/<PR#>.saturation.rst``
    - ``changes/<PR#>.skymatch.rst``
    - ``changes/<PR#>.source_catalog.rst``
    - ``changes/<PR#>.tweakreg.rst``
  </details>
- [ ] updated relevant tests
- [ ] updated relevant documentation
- [ ] updated relevant milestone(s)
- [ ] added relevant label(s)
- [ ] ran regression tests, post a link to the Jenkins job below. [How to run regression tests on a PR](https://github.com/spacetelescope/romancal/wiki/Running-Regression-Tests-Against-PR-Branches)
