import romancal.pipeline
import romancal.step
from romancal.stpipe import RomanPipeline, RomanStep
from romancal.stpipe.integration import get_steps


def test_get_steps():
    tuples = get_steps()

    assert {t[0].split(".")[-1] for t in tuples} == set(
        romancal.step.__all__ + romancal.pipeline.__all__
    )

    for class_name, class_alias, is_pipeline in tuples:
        # parse class path returning only class name
        _, name = class_name.rsplit(".", maxsplit=1)
        if is_pipeline:
            step_class = getattr(romancal.pipeline, name)
            assert issubclass(step_class, RomanPipeline)
        else:
            step_class = getattr(romancal.step, name)
        assert issubclass(step_class, RomanStep)
        assert step_class.class_alias == class_alias
