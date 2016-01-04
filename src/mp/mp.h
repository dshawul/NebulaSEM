#ifndef __MP_H
#define __MP_H

#include "mpi.h"
#include "my_types.h"

#if defined __DOUBLE
#	define MPI_SCALAR  MPI_DOUBLE
#else
#	define MPI_SCALAR  MPI_FLOAT
#endif

class MP {
public:
	enum {
		FIELD, END, FIELD_BLK
	};
	enum {
		OP_MAX, OP_MIN, OP_SUM, OP_PROD
	};
	MP(int argc,char* argv[]);
	~MP();
public:
	typedef MPI_Request REQUEST;

	static int n_hosts,host_id,name_len;
	static char host_name[PATH_MAX + 1];
	static int _start_time;
	static bool Terminated;
	static bool printOn;
	static char workingDir[PATH_MAX + 1];
	static void cleanup();
	static void loop();
	static void barrier();
	static int iprobe(int&,int&,int);
	static void send(int,int);
	static void recieve(int,int);
	static void printH(const char* format,...);
	static void print(const char* format,...);
	static bool hasElapsed(const Int);

	/*send and recieve messages*/
	template <class type>
	static void recieve(type* buffer,int size,int source,int message_id) {
		const int count = (size * sizeof(type) / sizeof(MPI_SCALAR));
		MPI_Recv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	template <class type>
	static void send(type* buffer,int size,int source,int message_id) {
		const int count = (size * sizeof(type) / sizeof(MPI_SCALAR));
		MPI_Send(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD);
	}
	template <class type>
	static void allreduce(type* sendbuf,type* recvbuf,int size, Int op) {
		const int count = (size * sizeof(type) / sizeof(MPI_SCALAR));
		MPI_Op mpi_op;
		switch(op) {
			case OP_MAX: mpi_op = MPI_MAX; break;
			case OP_MIN: mpi_op = MPI_MIN; break;
			case OP_SUM: mpi_op = MPI_SUM; break;
			case OP_PROD: mpi_op = MPI_PROD; break;
		}
		MPI_Allreduce(sendbuf,recvbuf,count,MPI_SCALAR,mpi_op,MPI_COMM_WORLD);
	}
	template <class type>
	static void irecieve(type* buffer,int size,int source,int message_id,void* request) {
		const int count = (size * sizeof(type) / sizeof(MPI_SCALAR));
		MPI_Irecv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
	}
	template <class type>
	static void isend(type* buffer,int size,int source,int message_id,void* request) {
		const int count = (size * sizeof(type) / sizeof(MPI_SCALAR));
		MPI_Isend(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
	}
	static void waitall(int count,void* request) {
		MPI_Waitall(count,(MPI_Request*)request,MPI_STATUS_IGNORE);
	}
};
#endif
