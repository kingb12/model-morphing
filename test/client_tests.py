import unittest
import Client
import Helpers
import Morph

morph = Helpers.make_morph()

class ClientTest(unittest.TestSuite):
    def test_translate_features():
        m2 = Client.translate_features(morph)
