#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sqlite3.h"

#include "spectral_database.h"


#define sql_check(err, db) { \
    if (err != SQLITE_OK) { \
        fprintf(stderr, "Error: %s\n", sqlite3_errmsg(db)); \
        return 1; \
    } \
}


int open_database(char const * path, sqlite3 ** connection)
{
    sqlite3 * db;
    sql_check(sqlite3_open(path, &db), db);
    if (db == NULL)
    {
        fprintf(stderr, "Error: failed to open %s.\n", path);
        return 1;
    }
    *connection = db;
    return 0;
}


int close_database(sqlite3 * connection)
{
    sql_check(sqlite3_close(connection), connection);
    return 0;
}


int compile_statement(sqlite3 * connection, char const * command, sqlite3_stmt ** statement)
{
    sqlite3_stmt * stmt;
    sql_check(sqlite3_prepare_v2(connection, command, -1, &stmt, NULL), connection);
    *statement = stmt;
    return 0;
}


int tips_data(sqlite3 * connection, int id, Tips_t * tips)
{
    /*Read tips data from the database.*/
    sqlite3_stmt * statement;
    char query[128];
    snprintf(query, 128,
             "select isotopologue_id, temperature, data from tips where molecule_id == %d",
             id);
    compile_statement(connection, query, &statement);
    int const n = 15*10000;
    tips->num_iso = 0;
    tips->temperature = malloc(sizeof(double)*n);
    tips->data = malloc(sizeof(double)*n);
    int i = 0;
    int current_iso = -1;
    while (sqlite3_step(statement) != SQLITE_DONE)
    {
        if (i >= n)
        {
            fprintf(stderr, "Error: buffer is too small, increase n.\n");
            return 1;
        }
        int iso = sqlite3_column_int(statement, 0);
        if (iso != current_iso)
        {
            tips->num_iso++;
            current_iso = iso;
        }
        tips->temperature[i] = sqlite3_column_double(statement, 1);
        tips->data[i] = sqlite3_column_double(statement, 2);
        i++;
    }
    if (current_iso == -1)
    {
       return -1;
    }
    tips->num_t = i / tips->num_iso;
    if (tips->num_t*tips->num_iso != i)
    {
        fprintf(stderr, "Error: tips data is not rectangular.\n");
        return 1;
    }
    sqlite3_finalize(statement);
    return 0;
}


/*Calculate total partition function.*/
double total_partition_function(Tips_t tips, double temperature, int iso)
{
    int i = iso*tips.num_t;
    double * t = tips.temperature + i; /*Pointer arithmetic.*/
    double * data = tips.data + i; /*Pointer arithemtic.*/
    i = (int)(floor(temperature)) - (int)(t[0]);
    return data[i] + (data[i+1] - data[i])*(temperature - t[i])/(t[i+1] - t[i]);
}


/*Reads the isotopologue masses.*/
int mass_data(sqlite3 * connection, int id, double * mass, int num_mass)
{
    sqlite3_stmt * statement;
    char query[128];
    snprintf(query, 128,
             "select isoid, mass from isotopologue where molecule_id == %d",
             id);
    compile_statement(connection, query, &statement);
    while (sqlite3_step(statement) != SQLITE_DONE)
    {
        int i = sqlite3_column_int(statement, 0);
        if (i == 0)
        {
            /*Weird HITRAN counting.*/
            i = 10;
        }
        if (i >= num_mass)
        {
            fprintf(stderr, "Error: buffer is too small, increase num_mass.\n");
            return 1;
        }
        mass[i-1] = sqlite3_column_double(statement, 1);
    }
    sqlite3_finalize(statement);
    return 0;
}


/*Reads the molecule id from the database.*/
int molecule_id(sqlite3 * connection, char * formula, int * id)
{
    /*Read tips data from the database.*/
    sqlite3_stmt * statement;
    char query[128];
    snprintf(query, 128,
             "select molecule from molecule_alias where alias == '%s'",
             formula);
    compile_statement(connection, query, &statement);
    *id = -1;
    while (sqlite3_step(statement) != SQLITE_DONE)
    {
        *id = sqlite3_column_int(statement, 0);
        break;
    }
    if (*id == -1)
    {
        fprintf(stderr, "Error: molecule %s not found in database.\n", formula);
        return 1;
    }
    sqlite3_finalize(statement);
    return 0;
}


/*Read the HITRAN parameters for a line.*/
int line_parameters(sqlite3_stmt * statement, LineParameter_t * parameter, double * mass)
{
    parameter->nu = sqlite3_column_double(statement, 0);
    parameter->sw = sqlite3_column_double(statement, 1);
    parameter->gamma_air = sqlite3_column_double(statement, 2);
    parameter->gamma_self = sqlite3_column_double(statement, 3);
    parameter->n_air = sqlite3_column_double(statement, 4);
    parameter->elower = sqlite3_column_double(statement, 5);
    parameter->delta_air = sqlite3_column_double(statement, 6);
    parameter->local_iso_id = sqlite3_column_int(statement, 7);
    if (parameter->local_iso_id == 0)
    {
        /*Weird HITRAN counting.*/
        parameter->local_iso_id = 10;
    }
    parameter->mass = mass[parameter->local_iso_id - 1];
    return 0;
}
