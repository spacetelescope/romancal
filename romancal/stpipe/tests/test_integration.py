from stpipe.utilities import import_class

from romancal.stpipe.integration import get_steps
from romancal.stpipe import RomanStep, RomanPipeline

import romancal.pipeline
import romancal.step


def test_get_steps():
    tuples = get_steps()

    assert {t[0].split(".")[-1] for t in tuples} == set(romancal.step.__all__ + romancal.pipeline.__all__)

    for class_name, class_alias, is_pipeline in tuples:
        step_class = import_class(class_name)
        assert issubclass(step_class, RomanStep)
        assert step_class.class_alias == class_alias
        if is_pipeline:
            assert issubclass(step_class, RomanPipeline)
