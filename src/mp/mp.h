#ifndef __MP_H
#define __MP_H

#include "mpi.h"
#include "tensor.h"

#if defined __DOUBLE
#	define MPI_SCALAR  MPI_DOUBLE
#else
#	define MPI_SCALAR  MPI_FLOAT
#endif

class MP {
public:
	enum {
		FIELD, END
	};
	MP(int argc,char* argv[]);
	~MP();
public:
	typedef MPI_Request REQUEST;

	static int n_hosts,host_id,name_len;
	static char host_name[512];
	static int _start_time;
	static bool Terminated;
	static void cleanup();
	static void loop();
	static void barrier();
	static int iprobe(int&,int&);
	static void send(int,int);
	static void recieve(int,int);
	static void printH(const char* format,...);
	static void print(const char* format,...);

	/*send and recieve messages*/
	template <class type>
	static void recieve(type* buffer,int size,int source,int message_id) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Recv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	template <class type>
	static void send(type* buffer,int size,int source,int message_id) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Send(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD);
	}
	template <class type>
	static void sendrecv(type* sbuffer,type* dbuffer,int size,int source,int dest, int message_id) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Status status;
		MPI_Sendrecv(sbuffer,count,MPI_SCALAR,dest  ,message_id,
					 dbuffer,count,MPI_SCALAR,source,message_id, 
					 MPI_COMM_WORLD,&status);
	}
	template <class type>
	static void allsum(type* sendbuf,type* recvbuf,int size) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Allreduce(sendbuf,recvbuf,count,MPI_SCALAR,MPI_SUM,MPI_COMM_WORLD);
	}
	template <class type>
	static void irecieve(type* buffer,int size,int source,int message_id,void* request) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Irecv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
	}
	template <class type>
	static void isend(type* buffer,int size,int source,int message_id,void* request) {
		const int count = (size * sizeof(type) / sizeof(Scalar));
		MPI_Isend(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
	}
	static void waitall(int count,void* request) {
		MPI_Waitall(count,(MPI_Request*)request,MPI_STATUS_IGNORE);
	}
};
#endif
