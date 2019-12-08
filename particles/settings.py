import os

BASE_DIR = os.path.dirname(__file__)
SAMPLE_DATA_FILE = os.path.join(BASE_DIR, 'data', 'sample.txt')


def _resource_path(name):
    return os.path.join(BASE_DIR, 'resources', name)


RECENTS_PATH = os.path.join(BASE_DIR, '.recents')
MAX_RECENTS = 10
MAX_RECENT_LENGTH = 45
