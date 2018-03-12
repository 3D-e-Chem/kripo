from kripodb.db import SqliteDb, SqliteDict


class LigandExpoDict(SqliteDict):
    def __init__(self, connection):
        super().__init__(connection, 'ligands', 'lig_id', 'mol')


class LigandExpoDb(SqliteDb):
    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS ligands (
            lig_id TEXT PRIMARY KEY,
            mol molblockgz
        )''')

    def as_dict(self):
        return LigandExpoDict(self.connection)
