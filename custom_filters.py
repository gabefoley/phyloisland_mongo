import phyloisland
from flask_admin.contrib.sqla.filters import BaseSQLAFilter


class GetUniqueSpecies(BaseSQLAFilter):
    def apply(self, query, value, alias="None"):

        species_list = []
        id_list = []

        for record in phyloisland.bio_db.values():
            species = (" ".join(record.annotations.get('organism').split()[0:2]))
            if species in species_list:
                continue
            else:
                species_list.append(species)
                id_list.append(record.id)
        # print(species_list)
        # print(id_list)
        return query.filter(self.get_column(alias).in_(id_list))

    def operation(self):
        return 'Yes'


class GetUniqueSpeciesSequenceRecord(BaseSQLAFilter):
    def apply(self, query, value, alias="None"):

        species_list = []
        id_list = []

        for record in phyloisland.bio_db.values():
            species = (" ".join(record.annotations.get('organism').split()[0:2]))
            if species in species_list:
                continue
            else:
                species_list.append(species)
                id_list.append(record.id)

        return query.filter(self.get_column(alias).in_(id_list))

    def operation(self):
        return 'Yes'