import argparse
import subprocess
import sys
from typing import Dict, List


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
    subprocess.run(command)


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
                "nphotons": 2000000000,
                "alpha": 5.,
                "n axicon": 1.45,
                "use_bottle": "true",
                "use_tracker": "false",
                "make_images": "false",
                "light_source": "point",
                "bottle_file": "clearBottle-large.params",
                "L2_file": "planoConvex.params",
                "L3_file": "achromaticDoublet.params",
                "image_source": "bessel-smear.dat",
                "data_folder": "settings"}

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
        if key == "wavelength":
            defaults[key] = str(defaults[key]).replace("e", "d")
        elif key in ["use_bottle", "use_tracker", "make_images"]:
            defaults[key] = str(defaults[key]).lower()

    # write out setting file
    with open("res/" + file, "w") as f:
        for key in defaults:
            f.write(str(defaults[key]) + " "*(35 - len(str(defaults[key]))) + "# " + key + "\n")
    return defaults


def create_spot_diags(bottles: List[List[str]]) -> None:
    """Summary

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
    """Summary

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


def offset_experiment() -> None:

    # default settings for image diagrams
    image_dict = {"light_source": "point", "use_tracker": "false", "make_images": "true",
                  "bottle_file": "clearBottle-large", "data_folder": "images-offset"}

    offsets = [i for i in range(4, 17, 2)]

    # run sims
    for i, offset in enumerate(offsets):
        image_dict["bottle_file"] = f"clearBottle-large_-{offset}mm.params"
        setup_f = f"test_{i}.params"
        make_settings(image_dict, setup_f)
        # run_sim(setup_f, image_dict)


def create_bessel_images(bottles: List[List[str]]) -> None:
    """Summary

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


parser = argparse.ArgumentParser(usage="%(prog)s [OPTION]", description="")
parser.add_argument("-s", "--spot", action="store_true", default=False,
                    help="Create spot diagrams.")
parser.add_argument("-p", "--point", action="store_true", default=False,
                    help="Create point/ring images.")
parser.add_argument("-b", "--bessel", action="store_true", default=False,
                    help="Create bessel/ring diagrams.")
parser.add_argument("-o", "--offset", action="store_true", default=False,
                    help="Run offset offset experiment on large bottle.")

bottles = [["clearBottle-large.params", "true"], ["clearBottle-small.params", "true"],
           ["clearBottle-ellipse.params", "true"], ["clearBottle-small.params", "false"]]

args = parser.parse_args()

if args.bessel:
    create_bessel_images(bottles)
if args.point:
    create_point_images(bottles)
if args.spot:
    create_spot_diags(bottles)
if args.offset:
    offset_experiment()
