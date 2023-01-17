#pragma once
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <windows.h>
#include <sqltypes.h>
#include <sql.h>
#include <sqlext.h>

namespace aqppp {
class SqlInterface {
 public:
  static int ConnectDb(SQLHANDLE &sqlconnectionhandle, std::string dsn,
                       std::string user, std::string pwd);
  static void ShowError(unsigned int handle_type, const SQLHANDLE &handle);
  static void ShowError(unsigned int handle_type, const SQLHANDLE &handle,
                        const std::string query_str);
  static void MakeSqlConnection(std::string odbc_name, std::string user_name,
                                std::string pwd,
                                SQLHANDLE &sqlconnectionhandle);
  /*
  return a query result of given string query.
  */
  static void SqlQuery(std::string query, SQLHANDLE &sqlstatementhandle);

  static std::string ComputeRandStr(SQLHANDLE &sqlConnectionHandle,
                                    std::string table_name, int seed,
                                    double sample_rate);

  // Create sample and small_sample table in MySQL database.
  static std::pair<double, double> CreateDbSamples(
      SQLHANDLE &sqlconnectionhandle, int seed, std::string db_name,
      std::string table_name, std::pair<double, double> sample_rates,
      std::pair<std::string, std::string> sample_names);
  static double CreateDbSample(SQLHANDLE &sqlconnectionhandle, int seed,
                               std::string db_name, std::string table_name,
                               double sample_rate, std::string sample_name);

  static double ReadDb(SQLHANDLE &sqlconnectionhandle,
                       std::vector<std::vector<double>> &o_table,
                       std::string db_name, std::string table_name,
                       std::string AGGREGATE_NAME,
                       std::vector<std::string> CONDITION_NAMES);
  static double Column2Numeric(SQLHANDLE &sqlstatementhandle, int col_id,
                               std::string col_name);
  static std::string ReadTableStr(std::string db_name, std::string table_name,
                                  std::string AGGREGATE_NAME,
                                  std::vector<std::string> CONDITION_NAMES);

  // Struct for storing names of numerical and categorical columns
  struct TableColumns {
    std::vector<std::string> all_columns;
    std::vector<std::string> numeric_columns;
    std::vector<std::string> categorical_columns;
  };

  // Returns TableColumns struct that lists the columns of the given table
  static const TableColumns GetTableColumns(SQLHANDLE &sqlConnectionHandle,
                                            std::string table_name);

 private:
};
}  // namespace aqppp
