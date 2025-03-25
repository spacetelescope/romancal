import time

from romancal.regtest.resource_tracker import ResourceTracker, TrackMemory, TrackRuntime


class FakeNode:
    def __init__(self):
        self.user_properties = []


class FakeRequest:
    def __init__(self):
        self.node = FakeNode()


def test_runtime():
    tracker = TrackRuntime()
    with tracker:
        time.sleep(1.0)
    assert abs(tracker.log()[1] - 1.0) < 0.01


def test_memory():
    tracker = TrackMemory()
    N = 1024 * 1024
    with tracker:
        b = b"0" * N  # noqa: F841
    assert abs(tracker.log()[1] - N) / N < 0.01


def test_resource_tracker():
    tracker = ResourceTracker()
    with tracker.track():
        pass
    fake_request = FakeRequest()
    tracker.log(fake_request)
    keys = {log[0] for log in fake_request.node.user_properties}
    assert keys == {"tracked-time", "tracked-peakmem"}


def test_log():
    tracker = ResourceTracker()
    fake_request = FakeRequest()
    with tracker.track(log=fake_request):
        pass
    keys = {log[0] for log in fake_request.node.user_properties}
    assert keys == {"tracked-time", "tracked-peakmem"}
