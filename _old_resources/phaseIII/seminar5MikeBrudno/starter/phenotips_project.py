import requests
from pprint import pprint as pp

###  Complete functions ###

def get_patients(url):
    """ Return the patient records (in JSON format) from the PhenoTips
    instance at url."""

    response = requests.get(url, headers={'Accept': 'application/json'})
    return response.json()


def get_patient(url, pid):
    """ Return the patient record (in JSON format) for patient pid from
    the PhentoTips instance at url."""

    response = requests.get(url + '/' + pid)
    return response.json()


def add_patient(url, patient):
    """ Add patient (in JSON format) to the PhenoTips instance at url. """
    response = requests.post(url, 
                            headers={'Content-Type': 'application/json'}, 
                            json=patient)
    return response


### Part 1: Your tasks ###

def get_patients_range(url, start, stop):
    """ Counting from the first existing record, return a JSON object with patients
    from the record start to the record stop for the PhenoTips instance at url.

    Note: some records from the PhenotIps instance may have been deleted, so the
    first and last records returned will likely have a higher IDs than the given
    start and stop positions.
    """

    # TODO


def delete_patient(url, pid):
    """ Delete patient with patient id pid from PhenotTips instance url.
    """

    # TODO


def get_phenotypic_info(patient):
    """ Given a patient's JSON representation, return a list of the phenotypic 
    features recorded for this patient.
    """

    # TODO


if __name__ == '__main__':

    url = 'https://playground.phenotips.org/rest/patients'

    # Get all patient records (response capped at 50 records) from the PhenoTips
    # playground in JSON format.
    # patients = get_patients(url)
    # pp(patients)

    # Print patient record for patient P0000055.
    # pp(get_patient(url, 'P0000065'))

    # Print the 10th to 15th records in the PhenoTips playground.
    # pp(get_patients_range(url, 10, 25))

    # patient_record = {"clinicalStatus": "affected",
    #  "patient_name": {"last_name": "Ng", "first_name": "Taylor"}, 
    #  "sex": "F", 
    #  "solved": {"status": "unsolved"} }
    # print(add_patient(url, patient_record))

    # Get the phenotypic features recorded for patient with pid 'P0000114'.
    # print(get_phenotypic_info(get_patient(url, 'P0000114')))

    # Delete a record (pick an id that exists in the PhenoTips playground), 
    # not 'P0006109' which is being used as a placeholder.
    # print(delete_patient(url, 'P0006109'))

