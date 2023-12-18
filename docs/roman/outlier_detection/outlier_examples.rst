Examples
========
Whether the data are contained in a list of ASDF files or provided as an ASN file, the
`ModelContainer` class must be used to properly handle the data that will be used in
the outlier detection step.

1. To run the outlier detection step (with the default parameters) on a list of 2 ASDF
   files named `"img_1.asdf"` and `"img_2.asdf"`:

        .. code-block:: python

                from romancal.outlier_detection import OutlierDetectionStep
                from romancal.datamodels import ModelContainer
                # read the file list into a ModelContainer object
                mc = ModelContainer(["img_1.asdf", "img_2.asdf"])
                step = OutlierDetectionStep()
                step.process(mc)

2. To run the outlier detection step (with the default parameters) on an ASN file
   called "asn_sample.json" with the following content:

        .. code-block:: json

                {
                        "asn_type": "None",
                        "asn_rule": "DMS_ELPP_Base",
                        "version_id": null,
                        "code_version": "0.9.1.dev28+ge987cc9.d20230106",
                        "degraded_status": "No known degraded exposures in association.",
                        "program": "noprogram",
                        "constraints": "No constraints",
                        "asn_id": "a3001",
                        "target": "none",
                        "asn_pool": "test_pool_name",
                        "products": [
                                {
                                "name": "files.asdf",
                                "members": [
                                        {
                                        "expname": "img_1.asdf",
                                        "exptype": "science"
                                        },
                                        {
                                        "expname": "img_2.asdf",
                                        "exptype": "science"
                                        }
                                ]
                                }
                        ]
                }

        .. code-block:: python

                from romancal.outlier_detection import OutlierDetectionStep
                from romancal.datamodels import ModelContainer
                # read the file list into a ModelContainer object
                mc = ModelContainer("asn_sample.json")
                step = OutlierDetectionStep()
                step.process(mc)

#. To run the outlier detection step (with the default parameters) on an ASN file
   called "asn_sample.json" (the files listed in the association file must have been
   processed through `TweakRegStep`) using the command line:

        .. code-block:: shell

                strun --disable-crds-steppar  --suffix='cal'  romancal.step.OutlierDetectionStep asn_sample.json
