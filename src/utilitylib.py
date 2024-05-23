#!/usr/bin/python3
from datetime import datetime, timezone  # ,timedelta
import yaml

# ---------------------------------------------------------------------------------------------------------------------


def load_yaml(yaml_file, yaml_path='.'):
    """
    Parameters:
    -----------
    yaml_file : name of yaml file
    """
    with open(f"{yaml_path}/{yaml_file}", 'r') as yhand:
        config_data = yaml.safe_load(yhand)
    return config_data

# ---------------------------------------------------------------------------------------------------------------------


def today_time_string():
    """
    Returns:
    -----------
    date string for plots
    """
    # Get the current UTC time
    current_time = datetime.now(timezone.utc)
    hour = current_time.hour  # Extract the hour in UTC

    # Determine the latest 6-hour interval
    if hour >= 18:
        utc_time = "18 UTC"
    elif hour >= 12:
        utc_time = "12 UTC"
    elif hour >= 6:
        utc_time = "06 UTC"
    else:
        utc_time = "00 UTC"

    # Return the formatted date and the appropriate 6-hour interval
    return current_time.strftime(f"%d.%m.%Y Time: {utc_time}")

# ---------------------------------------------------------------------------------------------------------------------


def station_number2strig(number_str):
    """
    Parameters:
    -----------
    number_str : WMO Station number as string
    """
    if number_str == '10035':
        station = 'Schleswig'
    elif number_str == '10113':
        station = 'Norderney'
    elif number_str == '10184':
        station = 'Greifswald'
    elif number_str == '10238':
        station = 'Bergen'
    elif number_str == '10304':
        station = '10304'
    elif number_str == '10393':
        station = 'Lindenberg'
    elif number_str == '10410':
        station = 'Essen'
    elif number_str == '10548':
        station = ' Meiningen'
    elif number_str == '10618':
        station = 'Idar-Oberstein'
    elif number_str == '10771':
        station = 'Kuemmersbruck'
    elif number_str == '10739':
        station = 'Stuttgart'
    elif number_str == '10868':
        station = 'Muenchen-Oberschlssheim'
    elif number_str == '06610':
        station = '06610'
    elif number_str == '11520':
        station = 'Praha-Libus'
    elif number_str == '11035':
        station = 'Wien'
    elif number_str == '16045':
        station = 'Rivolto'
    elif number_str == '11747':
        station = 'Prostejov'
    else:
        station = number_str
    return station
