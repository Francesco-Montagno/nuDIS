def pdg_id(part_name):
    mapping = {
        "d" : 1,   # down quark
        "u" : 2,   # up quark
        "s" : 3,   # strange quark
        "c" : 4,   # charm quark
        "b" : 5,   # bottom quark
        "t" : 6,   # top quark
        "db" : -1,   # down quark
        "ub" : -2,   # up quark
        "sb" : -3,   # strange quark
        "cb" : -4,   # charm quark
        "bb" : -5,   # bottom quark
        "tb" : -6,   # top quark
        "g" : 21,  # gluon
        "A" : 22,   # photon
        # Add more mappings as needed
    }
    return mapping.get(part_name, 999)  # Default to 999 if pdg_id is not in the mapping
