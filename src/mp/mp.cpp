#include <cstdarg>
#include "mp.h"
#include "system.h"

/*statics*/
int  MP::n_hosts;
int  MP::host_id;
int  MP::name_len;
char MP::host_name[512];
int  MP::_start_time = 0;
bool MP::Terminated = false;
bool MP::printOn = true;

/*Initialize*/
MP::MP(int argc,char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_hosts);
	MPI_Comm_rank(MPI_COMM_WORLD, &host_id);
	MPI_Get_processor_name(host_name, &name_len);
	_start_time = System::get_time();
	printf("Process [%d/%d] on %s : pid %d\n",
		host_id,n_hosts,host_name,System::get_pid());
	fflush(stdout);
}

/*finalize*/
MP::~MP() {
	MPI_Finalize();
}

/*send*/
void MP::send(int source,int message_id) {
	MPI_Send(MPI_BOTTOM,0,MPI_INT,source,message_id,MPI_COMM_WORLD);
}

/*recieve*/
void MP::recieve(int source,int message_id) {
	MPI_Recv(MPI_BOTTOM,0,MPI_INT,source,message_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

/*barrier*/
void MP::barrier() {
	MPI_Barrier(MPI_COMM_WORLD);
}

/*probe for messages*/
int MP::iprobe(int& source,int& message_id,int tag) {
	int flag;
	MPI_Status mpi_status;
	MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD,&flag,&mpi_status);
	if(flag) {
		message_id = mpi_status.MPI_TAG;
		source = mpi_status.MPI_SOURCE;
		return true;
	}
	return false;
}
/*print*/
void MP::printH(const char* format,...) {
	printf("%d [%d] ",System::get_time() - _start_time,host_id);
	va_list ap;
	va_start(ap, format);
	vprintf(format, ap);
	va_end(ap);
	fflush(stdout);
}
void MP::print(const char* format,...) {
	va_list ap;
	va_start(ap, format);
	vprintf(format, ap);
	va_end(ap);
	fflush(stdout);
}
/*cleanup*/
void MP::cleanup () {
	Terminated = true;
	printf("%d [%d] Exiting application\n", 
		System::get_time() - MP::_start_time, MP::host_id);
}
/*delay*/
bool MP::hasElapsed(const int delta) {
	static Int prev_time = 0;
	Int current_time;
	current_time = System::get_time();
	if(current_time - prev_time < delta)
		return false;
	prev_time = current_time;
	return true;
}