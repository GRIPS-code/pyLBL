from pkg_resources import get_entry_map, load_entry_point, Requirement
from re import match


plugins = get_entry_map("pyLBL")
models = plugins.keys()
molecular_lines = {key: value["Gas"].load() for key, value in plugins.items()}

cross_sections = {}
for key, value in plugins.items():
    if "CrossSection" in value:
        cross_sections[key] = value["CrossSection"].load()

collision_induced_absorption = {}
for key, value in plugins.items():
    if "CIA" in value:
        collision_induced_absorption[key] = value["CIA"].load()

continua = {}
for key, value in plugins.items():
    local_dict = {}
    for k, v in value.items():
        m = match(r"([A-Za-z0-9]+)Continuum", k)
        if m:
            local_dict[m.group(1)] = v.load()
    if local_dict:
        continua[key] = local_dict
