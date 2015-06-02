def postprocess_config(config):
    for key in "library species assembly targets sgrnas experiments samples".split():
        if key not in config:
            raise Exception("Missing key in config file: {}".format(key))
    for sample in config["samples"]:
        if not isinstance(config["samples"][sample], list):
            config["samples"][sample] = [config["samples"][sample]]
