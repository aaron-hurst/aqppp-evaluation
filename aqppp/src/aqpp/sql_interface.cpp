//#include"stdafx.h"
#include "tool.h"
#include "sql_interface.h"

namespace aqppp {

	void SqlInterface::ShowError(unsigned int handle_type, const SQLHANDLE& handle)
	{
		SQLWCHAR sqlstate[1024];
		SQLWCHAR message[1024];
		if (SQL_SUCCESS == SQLGetDiagRec(handle_type, handle, 1, sqlstate, NULL, message, 1024, NULL))
			std::wcout << "Message: " << message << "\nSQLSTATE: " << sqlstate << std::endl;
	}


	void SqlInterface::MakeSqlConnection(std::string odbc_name, std::string user_name, std::string pwd, SQLHANDLE &sqlconnectionhandle)
	{
		SQLHANDLE sqlenvhandle = NULL;
		if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &sqlenvhandle)) return;
		if (SQL_SUCCESS != SQLSetEnvAttr(sqlenvhandle, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0)) return;
		if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_DBC, sqlenvhandle, &sqlconnectionhandle)) return;

		switch (SQLConnect(sqlconnectionhandle, (SQLWCHAR*)odbc_name.c_str(), SQL_NTS, (SQLWCHAR*)user_name.c_str(), SQL_NTS, (SQLWCHAR*)pwd.c_str(), SQL_NTS))
		{
		case SQL_SUCCESS_WITH_INFO:
			std::cout << "success" << std::endl;
			ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
			break;
		default:
			std::cout << "error" << std::endl;
			ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
			getchar();
			return;
		}
	}

	int SqlInterface::ConnectDb(SQLHANDLE &sqlconnectionhandle, std::string dsn, std::string user, std::string pwd)
	{
		int retcode = 0;
		switch (SQLConnect(sqlconnectionhandle, (SQLWCHAR*)dsn.c_str(), SQL_NTS, (SQLWCHAR*)user.c_str(), SQL_NTS, (SQLWCHAR*)pwd.c_str(), SQL_NTS))
		{
		case SQL_SUCCESS_WITH_INFO:
			std::cout << "success" << std::endl;
			ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
			break;
		case SQL_INVALID_HANDLE:
		case SQL_ERROR:
			ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
			retcode = -1;
			break;
		default:
			ShowError(SQL_HANDLE_DBC, sqlconnectionhandle);
			break;
		}
		return retcode;
	}

	/*
	return a query result of given string query.
	*/
	void SqlInterface::SqlQuery(std::string query, SQLHANDLE &sqlstatementhandle)
	{
		std::wstring wquery = std::wstring(query.begin(), query.end());
		//SQLWCHAR wq[10000] = {};
		WCHAR* wq = const_cast<WCHAR*>(wquery.c_str());
		std::wcout << "sql_query: " << wquery << std::endl;
		if (SQL_SUCCESS != SQLExecDirect(sqlstatementhandle,wq, SQL_NTS))
		{
			ShowError(SQL_HANDLE_STMT, sqlstatementhandle);
		}
	}

	/*create sample and small_sample table in MySQL database.
	*/
	std::pair<double, double> SqlInterface::CreateDbSamples(SQLHANDLE &sqlconnectionhandle, int seed, std::string db_name, std::string table_name, std::pair<double, double> sample_rates, std::pair<std::string, std::string> sample_names)
	{
		double t1 = CreateDbSample(sqlconnectionhandle, seed+1, db_name, table_name, sample_rates.first, sample_names.first);
		double t2 = CreateDbSample(sqlconnectionhandle, seed+2, db_name, sample_names.first, sample_rates.second, sample_names.second);
		return { t1, t2 };
	}

	// Get column names for a given table and organise them into numerical and categorical
	SqlInterface::TableColumns SqlInterface::GetTableColumns(SQLHANDLE& sqlConnectionHandle, std::string table_name)
	{
		TableColumns results;
		SQLWCHAR buf_name[128];
		SQLWCHAR buf_type[128];
		std::vector<std::string> numeric_types{ "float", "int", "numeric" };  // TODO extend to include all numerical SQL types

		// Prepare SQL query to request all column names and types for the given table
		std::string stmt = "SELECT COLUMN_NAME, DATA_TYPE FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = '" + table_name + "';";
		std::wstring wstmt = std::wstring(stmt.begin(), stmt.end());
		const wchar_t* char_stmt = wstmt.c_str();
		SQLHANDLE h_stmt;
		RETCODE rc;
		rc = SQLAllocHandle(SQL_HANDLE_STMT, sqlConnectionHandle, &h_stmt);
		rc = SQLPrepare(h_stmt, (SQLWCHAR*)char_stmt, SQL_NTS);

		// Excecute the query
		rc = SQLExecDirect(h_stmt, (SQLWCHAR*)char_stmt, SQL_NTS);
		if (SQL_SUCCEEDED(rc))
		{
			SQLSMALLINT field_count = 0;
			SQLNumResultCols(h_stmt, &field_count);
			if (field_count == 2)
			{
				// Loop through the rows in the result set and parse into TableColumns results
				rc = SQLFetch(h_stmt);
				SQLLEN ret1;
				SQLLEN ret2;
				std::string column_name;
				std::string column_type;				
				while (SQL_SUCCEEDED(rc))
				{
					// Get data
					SQLGetData(h_stmt, 1, SQL_C_WCHAR, buf_name, sizeof(buf_name), &ret1);
					SQLGetData(h_stmt, 2, SQL_C_WCHAR, buf_type, sizeof(buf_type), &ret2);

					// Convert data to string
					if (ret1 <= 0) {
						column_name = std::string("(null");
					}
					else {
						std::wstring ws(buf_name);
						column_name = std::string(ws.begin(), ws.end());
					}
					if (ret2 <= 0) {
						column_type = std::string("(null");
					}
					else {
						std::wstring ws(buf_type);
						column_type = std::string(ws.begin(), ws.end());
					}

					// Push column name to correct vector
					if (std::find(std::begin(numeric_types), std::end(numeric_types), column_type) != std::end(numeric_types)) {
						results.numeric_columns.push_back(column_name);
					}
					else {
						results.categorical_columns.push_back(column_name);
					}

					// Get next row from query
					rc = SQLFetch(h_stmt);
				}
				rc = SQLFreeStmt(h_stmt, SQL_DROP);
			}
			else
			{
				std::cout << "Error: invalid number of fields returned: " << field_count << std::endl;
			}
		}
		else
		{
			std::cout << "error" << std::endl;
			ShowError(SQL_HANDLE_STMT, h_stmt);
		}
		return results;
	};

	std::string SqlInterface::ComputeRandStr(SQLHANDLE& sqlConnectionHandle, std::string table_name, int seed, double sample_rate)
	{
		TableColumns table_columns = GetTableColumns(sqlConnectionHandle, table_name);
		int tpseed = seed + 1;
		std::string numeric_part = std::to_string(seed);
		std::string catagory_part = "";
		for (std::string st : table_columns.numeric_columns)
		{
			numeric_part += ", " + std::to_string(tpseed) + "*" + st;
			tpseed++;
		}
		for (auto st : table_columns.categorical_columns)
		{
			catagory_part += ", " + st;
		}
		return "RAND(BINARY_CHECKSUM(" + numeric_part + catagory_part+"))<=" + std::to_string(sample_rate);
	}

	double SqlInterface::CreateDbSample(SQLHANDLE &sqlconnectionhandle, int seed, std::string db_name, std::string table_name, double sample_rate, std::string sample_name)
	{
		SQLHANDLE sqlstatementhandle = NULL;
		if (SQLAllocHandle(SQL_HANDLE_STMT, sqlconnectionhandle, &sqlstatementhandle) != SQL_SUCCESS) return -1;
		std::string sample_full_name = db_name + "." + sample_name;
		std::string table_full_name = db_name + "." + table_name;
		std::string drop_sample = "IF OBJECT_ID('" + sample_full_name + "', 'U') IS NOT NULL DROP TABLE " + sample_full_name + "; ";
		std::string create_sample = "SELECT * INTO " + sample_full_name + " FROM " + table_full_name + " WHERE " + ComputeRandStr(sqlconnectionhandle, table_name, seed, sample_rate) + ";";
		std::string create_sample_cstore_indx = "CREATE CLUSTERED COLUMNSTORE INDEX cci_"+sample_name+" ON " + sample_full_name + ";";
		std::cout << "Setting up sample database tables" << std::endl;
		std::cout << drop_sample << std::endl;
		std::cout << create_sample << std::endl;
		std::cout << create_sample_cstore_indx << std::endl;

		SqlQuery(drop_sample, sqlstatementhandle);
		double t1 = clock();
		SqlQuery(create_sample, sqlstatementhandle);
		SqlQuery(create_sample_cstore_indx, sqlstatementhandle);
		double create_sample_time = (clock() - t1) / CLOCKS_PER_SEC;
		SQLFreeHandle(SQL_HANDLE_STMT, sqlstatementhandle);
		return create_sample_time;
	}


	double SqlInterface::Column2Numeric(SQLHANDLE &sqlstatementhandle, int col_id, std::string col_name)
	{
		const int max_char_len = 300;
		transform(col_name.begin(), col_name.end(), col_name.begin(), [](unsigned char c) { return std::tolower(c); });
		if (col_name.find("date") != std::string::npos)  // date columns
		{
			double data = -1;
			TIMESTAMP_STRUCT ts;
			SQLGetData(sqlstatementhandle, col_id + 1, SQL_C_DATE, &ts, 0, NULL);
			data = ts.year * 10000 + ts.month * 100 + ts.day;
			return data;
		}
		if (col_name == "vendor_name")  // column only relevant for SSB dataset
		{
			double data = 0;
			char ts[10] = {};
			SQLGetData(sqlstatementhandle, col_id + 1, SQL_C_CHAR, ts, 10, NULL);
			std::string name = ts;
			if (name == "CMT") return 1;
			if (name == "DDS") return 2;
			if (name == "VTS") return 3;
			return -1;
		}
		if (col_name.find("time") != std::string::npos)  // time columns
		{
			double data = -1;
			TIMESTAMP_STRUCT ts;
			SQLGetData(sqlstatementhandle, col_id + 1, SQL_C_TYPE_TIMESTAMP, &ts, 0, NULL);
			data = ts.hour * 10000 + ts.minute * 100 + ts.second;
			return data;
		}
		
		// Otherwise: numerical data
		// NOTE: It was necessary to change the data type for the buffer in SQLGetData
		// to float due to the type of the data I am using.
		// TODO: Make this code adaptive so that it determine which type to use in SQLGetData
		// by itself (i.e. whether to use SQL_C_FLOAT, SQL_C_DOUBLE, SQL_C_INT, etc.).
		float data_float = 0;
		RETCODE rc;
		rc = SQLGetData(sqlstatementhandle, col_id + 1, SQL_C_FLOAT, &data_float, 0, NULL);
		if (!SQL_SUCCEEDED(rc))
		{
			std::cout << "Error getting data in SqlInterface::Column2Numeric." << std::endl;
			ShowError(SQL_HANDLE_STMT, sqlstatementhandle);
		}
		double data = (double)data_float;

		return data;
	}
	

	std::string SqlInterface::ReadTableStr(std::string db_name, std::string table_name, std::string AGGREGATE_NAME, std::vector<std::string> CONDITION_NAMES)
	{
		std::string demand_str = AGGREGATE_NAME;
		for (int i = 0; i < CONDITION_NAMES.size(); i++)
			demand_str += ", " + CONDITION_NAMES[i];

		std::string query = "SELECT " + demand_str + " FROM " + db_name + "." + table_name + ";";
		return query;
	}

	double SqlInterface::ReadDb(SQLHANDLE &sqlconnectionhandle, std::vector<std::vector<double>> &o_table, std::string db_name, std::string table_name, std::string AGGREGATE_NAME, std::vector<std::string> CONDITION_NAMES)
	{
		SQLHANDLE sqlstatementhandle = NULL;
		if (SQLAllocHandle(SQL_HANDLE_STMT, sqlconnectionhandle, &sqlstatementhandle) != SQL_SUCCESS) return -1;

		double st = clock();
		o_table = std::vector<std::vector<double>>();
		std::string query = ReadTableStr(db_name, table_name, AGGREGATE_NAME, CONDITION_NAMES);
		SqlInterface::SqlQuery(query, sqlstatementhandle);
		short int COL_NUM = 0;
		SQLNumResultCols(sqlstatementhandle, &COL_NUM);
		for (int i = 0; i < COL_NUM; i++)
			o_table.push_back(std::vector<double>());  // add a vector to the table for each column in the query
		while (SQLFetch(sqlstatementhandle) == SQL_SUCCESS)
		{
			double acc_data = Column2Numeric(sqlstatementhandle, 0, AGGREGATE_NAME);
			o_table[0].push_back(acc_data);
			for (int ci = 1; ci < COL_NUM; ci++)
			{
				double data = Column2Numeric(sqlstatementhandle, ci, CONDITION_NAMES[ci-1]);
				o_table[ci].push_back(data);
			}
		}
		double read_time = ((double)clock() - st) / CLOCKS_PER_SEC;
		SQLFreeHandle(SQL_HANDLE_STMT, sqlstatementhandle);
		return read_time;
	}
}