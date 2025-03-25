import time
import tracemalloc


class TrackRuntime:
    def __enter__(self):
        self._t0 = time.monotonic()

    def __exit__(self, exc_type, exc_value, traceback):
        self.value = time.monotonic() - self._t0

    def log(self):
        return ("tracked-time", self.value)


class TrackMemory:
    def __enter__(self):
        tracemalloc.start()

    def __exit__(self, exc_type, exc_value, traceback):
        _, self.value = tracemalloc.get_traced_memory()
        tracemalloc.stop()

    def log(self):
        return ("tracked-peakmem", self.value)


class ResourceTracker:
    def __init__(self):
        self.trackers = [TrackMemory(), TrackRuntime()]

    def __enter__(self):
        [t.__enter__() for t in self.trackers]

    def __exit__(self, exc_type, exc_value, traceback):
        [t.__exit__(exc_type, exc_value, traceback) for t in self.trackers]

    def log(self, user_properties):
        user_properties.extend(t.log() for t in self.trackers)


class ResourceTrackerManager:
    def __init__(self):
        self.named_trackers = {}

    def track(self, name):
        named_tracker = ResourceTracker()
        self.named_trackers[name] = named_tracker
        return named_tracker

    def log(self, user_properties, name=None):
        if name is None:
            self.named_trackers.popitem()[1].log(user_properties)
        else:
            self.named_trackers[name].log(user_properties)
