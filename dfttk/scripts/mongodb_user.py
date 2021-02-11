import random
import string

def get_random_string(length):
    letters = string.ascii_lowercase+string.ascii_uppercase+string.digits
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str

## Configuration

DB_HOST  = '146.186.149.69'
DB_PORT = 27018

## Imports
import click

passlength = 12

## Main function
def main():
    
    username = click.prompt("What is the new user's username?").strip()
    print()

    # Fireworks DB
    fws_db_name = username + '-fws'
    print(f"A database, {fws_db_name}, will be created with this user to store Fireworks.")
    fws_db_password = get_random_string(passlength)

    # Results DB, admin and read-only
    results_db_name = username + '-results'
    results_db_admin_username = username
    results_db_readonly_username = username + '-ro'
    print("\n")
    print(f"A database, {results_db_name}, will be created with this user to store results.")
    print("The database needs an admin and a read-only user.")
    results_db_admin_password = get_random_string(passlength)
    results_db_readonly_password= get_random_string(passlength)
    # Summary
    print("\n")
    print("The following actions will take place")
    print(f"Creating a Fireworks database, {fws_db_name}, with")
    print(f"    An admin user ({username}, password={fws_db_password})")
    print(f"Creating a Results database, {results_db_name}, with")
    print(f"    An admin user ({results_db_admin_username}, password={results_db_admin_password})")
    print(f"    A read-only user ({results_db_readonly_username}, password={results_db_readonly_password})")

    # Confirm
    """
    doit = click.confirm("Proceed with these actions?")
    if doit:
        print("Executing")
        raise NotImplementedError("Automatic user addition not implemented")
    else:
    """
    if True:
        print("Dry run")
        print("Run the following commands to make these changes:")
        print(f'use {fws_db_name}')
        print(f'db.createUser({{user: "{username}", pwd: "{fws_db_password}", roles: [{{role: "dbOwner", db: "{fws_db_name}"}}]}})')
        print(f'use {results_db_name}')
        print(f'db.createUser({{user: "{results_db_admin_username}", pwd: "{results_db_admin_password}", roles: [{{role: "dbOwner", db: "{results_db_name}"}}]}})')
        print(f'db.createUser({{user: "{results_db_readonly_username}", pwd: "{results_db_readonly_password}", roles: [{{role: "read", db: "{results_db_name}"}}]}})')

    print("Writing files")
    write_db_json(results_db_name, results_db_admin_username, results_db_admin_password, results_db_readonly_username, results_db_readonly_password, DB_HOST, DB_PORT)
    write_my_launchpad_yaml(fws_db_name, username, fws_db_password, DB_HOST, DB_PORT)


## Helper functions

TEMPLATE_DB_JSON = """{{
    "database": "{0}",
    "collection": "tasks",
    "admin_user": "{1}",
    "admin_password": "{2}",
    "readonly_user": "{3}",
    "readonly_password": "{4}",
    "host": "{5}",
    "port": {6},
    "aliases": {{}}
}}
"""

TEMPLATE_MY_LAUNCHPAD_YAML = """host: {0}
port: {1}
name: {2}
username: {3}
password: {4}
ssl_ca_file: null
strm_lvl: INFO
user_indices: []
wf_user_indices: []
"""

def write_db_json(dbname, adminuser, adminpwd, readonlyuser, readonlypwd, host, port, prefix=''):
    with open(prefix+'db.json', 'w') as fp:
        fp.write(TEMPLATE_DB_JSON.format(dbname, adminuser, adminpwd, readonlyuser, readonlypwd, host, port))

def write_my_launchpad_yaml(dbname, user, pwd, host, port, prefix=''):
    with open(prefix+'my_launchpad.yaml', 'w') as fp:
        fp.write(TEMPLATE_MY_LAUNCHPAD_YAML.format(host, port, dbname, user, pwd))

## Script
if __name__ == '__main__':
    main()


