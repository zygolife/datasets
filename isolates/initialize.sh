rm -f isolates.sqlite3
sqlite3 isolates.sqlite3 < schema.sql
sqlite3 isolates.sqlite3 < import.sql
