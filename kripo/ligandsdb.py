from kripodb.db import SqliteDb, SqliteDict


class LigandsDict(SqliteDict):
    def __init__(self, connection):
        super().__init__(connection, 'ligands', 'lig_id', 'mol')


class LigandsDb(SqliteDb):
    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS ligands (
            lig_id TEXT PRIMARY KEY,
            mol molblockgz
        )''')

    def as_dict(self):
        return LigandsDict(self.connection)
