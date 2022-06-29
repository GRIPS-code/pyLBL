#ifndef SPECTRAL_DATABASE_H_
#define SPECTRAL_DATABASE_H_

#include "sqlite3.h"


typedef struct LineParameter
{
    double nu;
    double sw;
    double gamma_air;
    double gamma_self;
    double n_air;
    double elower;
    double delta_air;
    int local_iso_id;
    double mass;
} LineParameter_t;


/*Parameters required to calculate the total parition function.*/
typedef struct Tips
{
    int num_iso;
    int num_t;
    double * temperature;
    double * data;
} Tips_t;


int open_database(char const * path, sqlite3 ** connection);


int close_database(sqlite3 * connection);


int compile_statement(sqlite3 * connection, char const * command,
                      sqlite3_stmt ** statement);


int tips_data(sqlite3 * connection, int id, Tips_t * tips);


/*Calculate total partition function.*/
double total_partition_function(Tips_t tips, double temperature, int iso);


/*Reads the isotopologue masses.*/
int mass_data(sqlite3 * connection, int id, double * mass, int num_mass);


/*Reads the molecule id from the database.*/
int molecule_id(sqlite3 * connection, char * formula, int * id);


/*Read the HITRAN parameters for a line.*/
int line_parameters(sqlite3_stmt * statement, LineParameter_t * parameter, double * mass);


#endif
