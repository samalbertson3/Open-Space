import sqlite3

def print_table(filename,tablename):
    c,cursor = connect(filename)
    for i in cursor.execute("""SELECT * FROM """+tablename):
        print i
    c.close()

def connect(filename):
    c = sqlite3.connect(filename)
    cursor = c.cursor()
    return c,cursor

def print_columns():
    c,cursor = connect('C:\Users\Sam\Desktop\orbits.db')
    li = []
    l = cursor.execute("""PRAGMA table_info(orbits);""")
    for i in l:
        li.append(i[1])
    return li
