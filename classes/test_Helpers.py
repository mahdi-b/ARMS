from Helpers import *
from nose.tools import assert_equals

class AttributeObject(object):
    def __init__(self, keys, values):
        for (key, value) in zip(keys, values):
            self.__dict__[key] = value

def test_mothur_buildOptionString():
    # Mothur option/attribute names (one filter and one update)
    filterName = "maxlength"
    updateName = "groups"
    # dummy values
    filterProp = 5
    updateProp = "a.group"
    # Test objects with various combinations of filter/update attributes
    emptyObj = AttributeObject([], [])
    filterObj = AttributeObject([filterName], [filterProp])
    updateObj = AttributeObject([updateName], [updateProp])
    completeObj = AttributeObject([filterName, updateName], [filterProp, updateProp])

    for object in [emptyObj, filterObj, updateObj, completeObj]:
        for filter in [True, False]:
            for update in [True, False]:
                try:
                    params = (object, filter, update)
                    attrDict = object.__dict__
                    # print(params)
                    testString = mothur_buildOptionString(*params)
                    print testString
                    # TODO actually check the testString for accuracy
                # if an update file was required, but the object didn't supply one
                except MissingMothurFileException:
                    if not updateName in attrDict:
                        assert True
                    else:
                        print("Unexpected MissingMothurFileException.")
                        assert False

                # if a filter was required but the object didn't supply one
                except MissingMothurFilterException:
                    if not filterName in attrDict:
                        assert True
                    else:
                        print("Unexpected MissingMothurFilterException.")
                        assert False

                # Something went wrong..
                except:
                    print sys.exc_info()
                    assert False

def test_buildCSL():
    testNames = ["a", "b", "c", "d"]
    x = AttributeObject(["a","b","c"],[0,1,None])
    xString = buildCSL(x, testNames)
    assert_equals(xString, "a=0, b=1")
    y = AttributeObject([],[])
    yString = buildCSL(y, testNames)
    assert_equals(yString, "")

test_buildCSL()
test_mothur_buildOptionString()