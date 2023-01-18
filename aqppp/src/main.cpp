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
  SQLHANDLE sqlenvhandle = NULL;
  SQLHANDLE sqlconnectionhandle = NULL;
  SQLHANDLE sqlstatementhandle = NULL;
  SQLWCHAR retconstring[1000];

  // Connect to SQL Server
  std::cout << "Attempting connection to SQL Server..." << std::endl;
  if (SQL_SUCCESS !=
      SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &sqlenvhandle))
    return -1;
  if (SQL_SUCCESS != SQLSetEnvAttr(sqlenvhandle, SQL_ATTR_ODBC_VERSION,
                                   (SQLPOINTER)SQL_OV_ODBC3, 0))
    return -1;
  if (SQL_SUCCESS !=
      SQLAllocHandle(SQL_HANDLE_DBC, sqlenvhandle, &sqlconnectionhandle))
    return -1;
  switch (SQLDriverConnect(
      sqlconnectionhandle, NULL,
      (SQLWCHAR*)L"DRIVER={SQL "
                 L"Server};SERVER=localhost,1434;DATABASE=uci_household_power_"
                 L"consumption;Trusted=true;",
      SQL_NTS, retconstring, 1024, NULL, SQL_DRIVER_NOPROMPT)) {
    case SQL_SUCCESS_WITH_INFO:
      std::cout << "success" << std::endl;
      aqppp::SqlInterface::ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
      break;
    default:
      std::cout << "error" << std::endl;
      aqppp::SqlInterface::ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
      return -1;
  }

  // Run experiments
  if (RUN_DEMO) expDemo::ComprehensiveExp::Exp(sqlconnectionhandle);
  if (RUN_COMPARISON)
    exp_comparison::ComparisonExperiment(sqlconnectionhandle).RunExperiment();

  std::cout << "All experiments completed." << std::endl;
  SQLDisconnect(sqlconnectionhandle);
  SQLFreeHandle(SQL_HANDLE_DBC, sqlconnectionhandle);
  SQLFreeHandle(SQL_HANDLE_ENV, sqlenvhandle);
  return 0;
}
