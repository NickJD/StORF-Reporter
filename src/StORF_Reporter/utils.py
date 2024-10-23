import collections

def sortORFs(tool_ORFs):  # Will only sort by given start position
    tool_ORFs_Sorted = sorted(tool_ORFs.items(), key=lambda v: int(v[0].split(",")[0]))
    tool_ORFs_Sorted = collections.OrderedDict(tool_ORFs_Sorted)
    return tool_ORFs_Sorted