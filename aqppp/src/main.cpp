#define NOMINMAX
#define _NO_CRT_STDIO_INLINE
#pragma comment(linker, "/STACK:3000000000")  // for batch_read in
                                              // precomputation
#pragma comment(linker, "/HEAP:3000000000")  // for batch_read in precomputation

#include <Windows.h>
#include <sqlext.h>
#include <stdio.h>

#include <any>
#include <cstdio>
#include <exception>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>

#include "aqpp\configuration.h"
#include "aqpp\sql_interface.h"
#include "exp\comparison.h"
#include "exp\comprehensive_exp.h"

#define RUN_DEMO 0
#define RUN_COMPARISON 1

int main() {
  // TODO new instructions
  // srand(0);
  // std::cout <<
  // "---------------------------------------------------------------------------------------------------"
  // << std::endl; std::cout << "The parameters are set at structure 'Settings'
  // in comment_content.h. This code was written using visual studio 2017, and
  // run in Windows System." << std::endl; std::cout << "It uses ODBC to connect
  // SQL Server 2017. The connection string could be found in the main function
  // (you should set or create related account and ODBC). Please change the
  // connection string based on your setting. " << std::endl; std::cout <<
  // "Besides, please load lineitem table into SQLServer, using the name in
  // parameter. For example, if you set par.DB_NAME = 'skew_s100_z2.dbo', and
  // par.TABLE_NAME='lineitem', then create a database called skew_s100_z2, and
  // a table called lineitem in SQLServer. " << std::endl; std::cout<<"Then load
  // lineitem.tbl (generated using TPCD benchmark by command -s 100 -z 2 -T L)
  // into the lineitem table in SQLServer. Then the code should work. The
  // results will be saved on exp_result folder." << std::endl; std::cout <<
  // "---------------------------------------------------------------------------------------------------"
  // << std::endl; std::cout << "To run the demo, input 1" << std::endl;
  ////std::cout << "2--query generation" << std::endl;

  // Define SQL handles
  SQLHANDLE sql_env_handle = NULL;
  SQLHANDLE sql_connection_handle = NULL;
  SQLWCHAR retconstring[1000];

  // Connect to SQL Server
  std::cout << "Attempting connection to SQL Server..." << std::endl;
  if (SQL_SUCCESS !=
      SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &sql_env_handle))
    return -1;
  if (SQL_SUCCESS != SQLSetEnvAttr(sql_env_handle, SQL_ATTR_ODBC_VERSION,
                                   (SQLPOINTER)SQL_OV_ODBC3, 0))
    return -1;
  if (SQL_SUCCESS !=
      SQLAllocHandle(SQL_HANDLE_DBC, sql_env_handle, &sql_connection_handle))
    return -1;
  switch (SQLDriverConnect(
      sql_connection_handle, NULL,
      (SQLWCHAR*)L"DRIVER={SQL Server};"
                 L"SERVER=localhost,1434;"
                 L"DATABASE=uci_household_power_consumption;"
                 L"Trusted=true;",
      SQL_NTS, retconstring, 1024, NULL, SQL_DRIVER_NOPROMPT)) {
    case SQL_SUCCESS_WITH_INFO:
      std::cout << "Success!" << std::endl;
      aqppp::SqlInterface::ShowError(SQL_HANDLE_DBC, sql_connection_handle);
      break;
    default:
      std::cout << "Error." << std::endl;
      aqppp::SqlInterface::ShowError(SQL_HANDLE_DBC, sql_connection_handle);
      return -1;
  }

  // Run experiments
  std::cout << "Beginning experiments..." << std::endl;
  if (RUN_DEMO) expDemo::ComprehensiveExp::Exp(sql_connection_handle);
  if (RUN_COMPARISON)
    exp_comparison::ComparisonExperiment(sql_connection_handle).RunExperiment();

  std::cout << "All experiments completed." << std::endl;
  SQLDisconnect(sql_connection_handle);
  SQLFreeHandle(SQL_HANDLE_DBC, sql_connection_handle);
  SQLFreeHandle(SQL_HANDLE_ENV, sql_env_handle);
  return 0;
}
