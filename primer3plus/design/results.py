import re


def parse_primer3_results(results_dict):
    """Parse the primer3 results.

    :param results_dict: :type results_dict: :return: :rtype:
    """
    num_pairs = results_dict["PRIMER_PAIR_NUM_RETURNED"]
    num_left = results_dict["PRIMER_LEFT_NUM_RETURNED"]
    num_right = results_dict["PRIMER_RIGHT_NUM_RETURNED"]

    pairs = {}
    other = {}
    for i in range(max([num_pairs, num_left, num_right])):
        pairs.setdefault(i, {})
    key_pattern = r"PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)_(?P<key>.+)"
    location_pattern = r"PRIMER_(?P<label>[a-zA-Z]+)_(?P<pair_id>\d+)\s*$"
    for k in results_dict:
        m = re.match(key_pattern, k)
        loc_m = re.match(location_pattern, k)
        if m:
            groupdict = m.groupdict()
            pair_id = int(groupdict["pair_id"])
            label = groupdict["label"]
            key = groupdict["key"]
            pairdict = pairs[pair_id]
            pairdict.setdefault(label, {})
            pairdict[label][key] = results_dict[k]
        elif loc_m:
            groupdict = loc_m.groupdict()
            pair_id = int(groupdict["pair_id"])
            label = groupdict["label"]
            pairdict = pairs[pair_id]
            pairdict.setdefault(label, {})
            pairdict[label]["location"] = results_dict[k]
        else:
            other[k] = results_dict[k]
    return pairs, other
