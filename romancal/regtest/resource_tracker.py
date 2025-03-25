"""Resource tracking regtest utilities

Can be used within module-scoped fixtures (often used to
run Steps or Pipelines) or within tests.

For uses where the resource usage occurs within a test:

>>>
>> def test_long_step(resource_tracker, request):
>>     with resource_tracker.track(log=request):
>>         # something that takes memory and time
>>         pass

For a module-scoped fixture the resource tracking can
be performed in the fixture but the logging/reporting of
the resource usage must occur during a test:

>>>
>> @pytest.fixture(scope="module")
>> def resource_tracker():
>>     return ResourceTracker()
>>
>> @pytest.fixture()
>> def log_tracked_resources(resource_tracker, request):
>>     def callback():
>>         resource_tracker.log(request)
>>
>>     yield callback
>>
>> @pytest.fixture
>> def my_long_fixture(resource_tracker):
>>     with resource_tracker.track():
>>         # something that takes memory and time
>>         pass
>>
>> def test_log_tracked_resources(log_tracked_resources, my_long_fixture):
>>     log_tracked_resources()

Use of the module-scoped fixture has fixture-reuse
considerations similar to the ``rtdata_module`` fixture. Having
more than one module scoped fixture that uses ``resource_tracker``
per module is discouraged (as both will use the same ``ResourceTracker``
instance). Parameterization of a fixture using ``resource_tracker``
is supported (same as ``rtdata_module``).
"""

import time
import tracemalloc
from contextlib import ExitStack, contextmanager


class TrackRuntime:
    """Runtime tracker context."""

    def __enter__(self):
        self._t0 = time.monotonic()

    def __exit__(self, exc_type, exc_value, traceback):
        self.value = time.monotonic() - self._t0

    def log(self):
        return ("tracked-time", self.value)


class TrackPeakMemory:
    """Peak memory tracker context."""

    def __enter__(self):
        tracemalloc.start()

    def __exit__(self, exc_type, exc_value, traceback):
        _, self.value = tracemalloc.get_traced_memory()
        tracemalloc.stop()

    def log(self):
        return ("tracked-peakmem", self.value)


class ResourceTracker:
    """Track resources used during track context."""

    def __init__(self):
        self._trackers = [TrackPeakMemory(), TrackRuntime()]

    def log(self, request):
        """Log tracked resource usage to the pytest request user properties.

        Parameters
        ----------
        request : pytest.FixtureRequest
            Must be a function-scoped pytest request fixture result.
        """
        request.node.user_properties.extend(t.log() for t in self._trackers)

    @contextmanager
    def track(self, log=None):
        """Context during which resources are tracked.

        Parameters
        ----------
        log : pytest.FixtureRequest, optional
            If provided, log the usage to the provided request fixture result.
        """
        try:
            with ExitStack() as stack:
                [stack.enter_context(t) for t in self._trackers]
                yield self
        finally:
            if log:
                self.log(log)
