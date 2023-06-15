(* ::Package:: *)

Needs["DatabaseLink`"];
JDBCDrivers["MySQL(Connector/J)"];

V27SQLConn[]:=Module[{},
	OpenSQLConnection[JDBC["MySQL(Connector/J)","veadbs-03:3306/v_27"],"Username"->"ebra","Password"->"eb125ra"]
]
HovsoreSQLConn[]:=Module[{},
	OpenSQLConnection[JDBC["MySQL(Connector/J)","veadbs-03:3306/hovsore"],"Username"->"ebra","Password"->"eb125ra"]
]
HovsoreTSSQLConn[]:=Module[{},
	OpenSQLConnection[JDBC["MySQL(Connector/J)","veadbs-03:3307/hovsore_20hz"],"Username"->"ebra","Password"->"eb125ra"]
]
