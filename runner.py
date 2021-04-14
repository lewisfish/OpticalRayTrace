import argparse
import subprocess
import sys
from typing import Dict, List

import numpy as np


def save_settings(file: str, place: str):
    """Save a record of the settings file to the data folder of
       current experiment.

    Parameters
    ----------
    file : str
        settings file to copy and rename.
    place : str
        Place to copy settings file to.
    """

    command = f"cp {file} data/{place}/settings.params"
    command = command.split()
    subprocess.run(command)


def run_sim(settings_file: str, settings_dict: Dict) -> None:
    """Summary

    Parameters
    ----------
    settings_file : str
        Settings file to pass to simulation
    settings_dict : Dict
        settings dictionary

    """

    # set ncores
    if settings_dict["use_tracker"] == "true":
        ncores = 1
    else:
        ncores = 32

    # run simulatrion
    command = f"./install.sh -n {ncores} -f {settings_file}"
    command = command.split(" ")
    subprocess.run(command, check=True)


def make_settings(user_dict: Dict, file: str) -> None:
    """Create file continang the current simulations settings

    Parameters
    ----------
    user_dict : Dict
        User defined keys which overide the default values
    file : str
        setting file to save parameters in

    Returns
    -------
    None

    """

    # default values
    defaults = {"ring_width": 0.5e-3,
                "wavelength": "785d-9",
                "nphotons": 1000000000,
                "alpha": 5.,
                "n axicon": 1.45,
                "use_bottle": "true",
                "use_tracker": "false",
                "make_images": "false",
                "image_diameter": 1e-2,
                "fibre_offset": 0.0,
                "light_source": "point",
                "iris": "none",
                "iris_size": 1.0,  # factor of lens radius
                "bottle_file": "clearBottle-large.params",
                "L2_file": "planoConvex-f39.9mm.params",
                "L3_file": "achromaticDoublet-f50.0mm.params",
                "image_source": "bessel-smear.dat",
                "data_folder": "settings",
                "isors_offset": 0.0}

    # check all user keys are valid
    for user_key in user_dict:
        if user_key in defaults:
            defaults[user_key] = user_dict[user_key]
        else:
            print(f"No such key: {user_key} with value {user_dict[user_key]}")
            sys.exit(0)

    # need to process user values:
    # convert e to d
    # lower case booleans
    for key in defaults:
        if key == "wavelength" or key == "image_diameter":
            defaults[key] = str(defaults[key]).replace("e", "d")
        elif key in ["use_bottle", "use_tracker", "make_images", "iris"]:
            defaults[key] = str(defaults[key]).lower()

    # write out setting file
    with open("res/" + file, "w") as f:
        for key in defaults:
            f.write(str(defaults[key]) + " "*(35 - len(str(defaults[key]))) + "# " + key + "\n")

    save_settings("res/" + file, defaults["data_folder"])

    return defaults


def create_spot_diags(bottles: List[List[str]]) -> None:
    """ -s setting

    Parameters
    ----------
    bottles : List[List[str]]
        list of settings files for which bottle to use in sim.
    """

    # default settings for spot diagrams
    spot_dict = {"nphotons": 100, "use_tracker": "true", "light_source": "spot",
                 "bottle_file": "clearBottle-small.params",
                 "data_folder": "spot-diag"}

    # run sims
    for i, setting in enumerate(bottles):
        spot_dict["bottle_file"] = setting[0]
        spot_dict["use_bottle"] = setting[1]
        setup_f = f"test_{i}.params"
        make_settings(spot_dict, setup_f)
        run_sim(setup_f, spot_dict)


def create_point_images(bottles: List[List[str]]) -> None:
    """ -p setting

    Parameters
    ----------
    bottles : List[List[str]]
        list of settings files for which bottle to use in sim.
    """

    # default settings for image diagrams
    image_dict = {"light_source": "point", "use_tracker": "false", "make_images": "true",
                  "bottle_file": "clearBottle-small.params", "data_folder": "images"}

    # run sims
    for i, setting in enumerate(bottles):
        image_dict["bottle_file"] = setting[0]
        image_dict["use_bottle"] = setting[1]
        setup_f = f"test_{i}.params"
        make_settings(image_dict, setup_f)
        run_sim(setup_f, image_dict)


def iris_experiment(bottles: List[List[str]]) -> None:
    """ -i setting

    Parameters
    ----------
    bottles : List[List[str]]
        list of settings files for which bottle to use in sim.
    """

    # default settings for image diagrams
    image_dict = {"light_source": "point", "use_tracker": "false", "make_images": "true",
                  "bottle_file": "clearBottle-small.params", "data_folder": "iris"}

    irises = ["before", "after", "none"]
    sizes = [1.0, 0.8, 0.6, 0.4, 0.2]
    # run sims
    for i, setting in enumerate(bottles):
        for iris in irises:
            for size in sizes:
                image_dict["bottle_file"] = setting[0]
                image_dict["use_bottle"] = setting[1]
                image_dict["iris"] = iris
                image_dict["iris_size"] = size
                setup_f = f"test_{i}_{iris}_{size}.params"
                make_settings(image_dict, setup_f)
                run_sim(setup_f, image_dict)
                if iris == "none":
                    # only need to do one size of iris for no iris as it dosent make any sense really...
                    break


def offset_experiment() -> None:
    """-o setting

    """

    # default settings for image diagrams
    image_dict = {"light_source": "point", "use_tracker": "false",
                  "make_images": "true", "bottle_file": "clearBottle-large",
                  "data_folder": "images-offset"}

    offsets = [i for i in range(4, 17, 2)]

    # run sims
    for i, offset in enumerate(offsets):
        image_dict["bottle_file"] = f"clearBottle-large_-{offset}mm.params"
        setup_f = f"test_{i}.params"
        make_settings(image_dict, setup_f)
        run_sim(setup_f, image_dict)


def create_bessel_images(bottles: List[List[str]]) -> None:
    """-b option

    Parameters
    ----------
    bottles : List[List[str]]
        list of settings files for which bottle to use in sim.
    """

    # default settings for image diagrams
    image_dict = {"light_source": "image", "use_tracker": "false", "make_images": "true",
                  "bottle_file": "clearBottle-small.params", "data_folder": "images"}

    # run sims
    for i, setting in enumerate(bottles):
        image_dict["bottle_file"] = setting[0]
        image_dict["use_bottle"] = setting[1]
        setup_f = f"test_{i}.params"
        make_settings(image_dict, setup_f)
        run_sim(setup_f, image_dict)


def lens_experiment(bottles: List[List[str]]) -> None:
    """-l setting

    Parameters
    ----------
    bottles : List[List[str]]
        Description

    """

    # default settings for image diagrams
    image_dict = {"light_source": "point", "use_tracker": "false",
                  "make_images": "false", "bottle_file": "clearBottle-large",
                  "data_folder": "images-lens"}

    # L2 focal lengths
    L2_flengths = ["59.8", "49.8", "39.9", "34.9", "29.9"]
    # L3 focal lengths
    L3_flengths = ["40.0", "45.0", "50.0", "60.0", "75.0"]

    # run sims
    for k, L3focal in enumerate(L3_flengths):
        image_dict["L3_file"] = f"achromaticDoublet-f{L3focal}mm.params"
        for j, L2focal in enumerate(L2_flengths):
            image_dict["L2_file"] = f"planoConvex-f{L2focal}mm.params"
            for i, bottle in enumerate(bottles):
                image_dict["bottle_file"] = bottle[0]
                image_dict["use_bottle"] = bottle[1]
                setup_f = f"test_{i}_{j}_{k}.params"
                make_settings(image_dict, setup_f)
                run_sim(setup_f, image_dict)


def bessel_params_experiment() -> None:
    # TODO once bpm.py is fixed
    pass
def iSORS_vs_Bessel():
    """isb option"""

    param_dict = {"light_source": "point", "use_tracker": "false",
                  "make_images": "true", "image_diameter": "1e-2",
                  "data_folder": "iSORS_vs_Bessel", "use_bottle": "true",
                  "ring_width": 0.5e-3, "alpha": 5, "n axicon": 1.45,
                  "L2_file": "planoConvex-f39.9mm.params"}

    # experiment parameters
    sources = ["isors", "point"]
    offsets = np.linspace(0., 1.5e-3, 7)

    # parameters required to calculate bottle offset for bessel beam
    n_axicon = param_dict["n axicon"]
    alpha = param_dict["alpha"] * np.pi / 180.
    file = param_dict["L2_file"]
    # read lens 2 file to get its back focal distance
    with open("res/" + file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "fb" in line:
                data = line
                break
    items = data.split(" ")
    l2fb = float(items[0].replace("d", "e"))
    init_dist = 97.3e-3  # distance from axicon to L1

    # read bottle file to get radius
    with open("res/" + "clearBottle-small_0.0mm.params") as f:
        lines = f.readlines()
        for line in lines:
            if "radius a" in line:
                data = line
                break
    items = data.split(" ")
    bottle_radiusA = float(items[0].replace("d", "e"))

    # carry out experiment
    for source in sources:
        param_dict["light_source"] = source
        for offset in offsets:
            param_dict["isors_offset"] = offset
            if source == "isors":
                param_dict["bottle_file"] = f"clearBottle-small_0.0mm.params"
            else:
                # calculate distance to move bottle to get required iSORS spatial offset
                prop_offset = ((l2fb*(offset + param_dict["ring_width"])) / (init_dist*np.tan(alpha* (n_axicon - 1)))) - bottle_radiusA
                create_bottle_file(prop_offset)
                param_dict["bottle_file"] = f"clearBottle-small_iSORS.params"

            setup_f = f"test.params"
            make_settings(param_dict, setup_f)
            run_sim(setup_f, param_dict)


def create_bottle_file(offset: float):
    """
    Create bottle file for iSORS experiment comparison with required offset of bottle to match iSORS setup

    Parameters
    ----------
    offset : float
        Required offset.
    """

    lines = ["2.d-3           thickness mm",
             "17.5d-3         radius a mm in z direction",
             "17.5d-3         radius b mm in y direction. if radius a != radius b then bottle is elliptical",
             "0.0             x",
             "0.0             y bottle position, centre of",
             f"{offset}             z",
             "1.5130          b1 glass dispersion coiefficents (soda-lime clear) https://refractiveindex.info/?shelf=glass&book=soda-lime&page=Rubin-clear",
             "0.003169        b2",
             "0.003962        b3",
             "1.35265         c1 bottle contents coifficents (alcohol) https://refractiveindex.info/?shelf=organic&book=ethanol&page=Rheims",
             "0.00306         c2",
             "0.00002         c3"]

    with open("res/clearBottle-small_iSORS.params", "w") as f:
        for line in lines:
            f.write(line + "\n")


parser = argparse.ArgumentParser(usage="%(prog)s [OPTION]", description="")
parser.add_argument("-s", "--spot", action="store_true", default=False,
                    help="Create spot diagrams.")
parser.add_argument("-p", "--point", action="store_true", default=False,
                    help="Create point/ring images.")
parser.add_argument("-b", "--bessel", action="store_true", default=False,
                    help="Create bessel/ring diagrams.")
parser.add_argument("-o", "--offset", action="store_true", default=False,
                    help="Run offset experiment on large bottle.")
parser.add_argument("-i", "--iris", action="store_true", default=False,
                    help="Run iris experiment on bottles.")
parser.add_argument("-l", "--lens", action="store_true", default=False,
                    help="Run lens experiments.")
parser.add_argument("-bp", "--bessel_params", action="store_true", default=False,
                    help="Run bessel parameters experiment.")

parser.add_argument("-a", "--all", action="store_true", default=False,
                    help="Run all experiments.")


bottles = [["clearBottle-large.params", "true"], ["clearBottle-small.params", "true"],
           ["clearBottle-ellipse.params", "true"], ["clearBottle-small.params", "false"]]

args = parser.parse_args()

if args.bessel or args.all:
    create_bessel_images(bottles)
if args.point or args.all:
    create_point_images(bottles)
if args.spot or args.all:
    create_spot_diags(bottles)
if args.offset or args.all:
    offset_experiment()
if args.iris or args.all:
    iris_experiment(bottles)
if args.lens or args.all:
    bottles = [["clearBottle-large.params", "true"], ["clearBottle-small.params", "true"],
               ["clearBottle-ellipse.params", "true"]]
    lens_experiment(bottles)
