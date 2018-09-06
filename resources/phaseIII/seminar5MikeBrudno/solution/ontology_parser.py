class Ontology:
    
    def __init__(self, ontology_file):
        """ Initialize a new Ontology object with data from ontology_file."""

        # Mapping of phenotypic ids to names.
        self.id_to_name = {}

        # Mapping of phenotypic ids to parents in the
        # human phenotype ontology.
        self.id_to_parents = {}

        # human_phenotype onology file to populate
        # the dictionaries id_to_name and id_to_parents.
        self.populate_dictionaries(ontology_file)

    def populate_dictionaries(self, ontology_file):
        """ Populate the id_to_name and id_to_parents dictionaries with
        data from ontology_file."""

        for line in ontology_file:
            if line.startswith('id:'):
                id_num = line.strip().split(':')[-1].strip()
                line = ontology_file.readline()
                name = ''
                parents = []
                while line != '\n':
                    if line.startswith('name:'):
                        name = line.strip().split(':')[-1].strip()
                    elif line.startswith('is_a:'):
                        parent = line.strip().split(':')[-1].split(' ! ')[0].strip()
                        parents.append(parent)
                    line = ontology_file.readline()
                self.id_to_parents[id_num] = parents
                self.id_to_name[id_num] = name


    def get_id_to_name(self):
        """ Return the id_to_name dictionary. """

        return self.id_to_name


    def get_id_to_parents(self):
        """ Return the id_to_parents dictionary."""

        return self.id_to_parents


    def get_parents(self, id_num):
        """ Return a list of the parents of id_num. """

        return self.id_to_parents[id_num]


    def get_name(self, id_num):
        """ Return the name associated with id_num."""

        return self.id_to_name[id_num]