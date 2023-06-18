#ifndef __MP_H
#define __MP_H

#include "mpi.h"
#include "types.h"

#if defined USE_DOUBLE
#   define MPI_SCALAR  MPI_DOUBLE
#else
#   define MPI_SCALAR  MPI_FLOAT
#endif

#define  MY_PATH_MAX  512

/**
  Class for multi-processor support via MPI
 */
class MP {
    public:
        /** Message types */
        enum {
            FIELD,      /**< Field data marker */
            END,        /**< END of communicatins marker */
            FIELD_BLK   /**< Field data marker when not in iteration*/
        };
        /** Global reduction types */
        enum {
            OP_MAX, /**< Global maximum */
            OP_MIN, /**< Global minimum */
            OP_SUM, /**< Global sum */
            OP_PROD /**< Global product */
        };
        MP(int argc,char* argv[]);
        ~MP();
    public:
        typedef MPI_Request REQUEST;

        static int n_hosts,host_id,name_len;
        static char host_name[MY_PATH_MAX + 1];
        static int _start_time;
        static bool Terminated;
        static bool printOn;
        static char workingDir[MY_PATH_MAX + 1];
        static void cleanup();
        static void loop();
        static void barrier();
        static int iprobe(int&,int&,int);
        static void send(int,int);
        static void recieve(int,int);
        static void printH(const char* format,...);
        static void print(const char* format,...);
        static bool hasElapsed(const Int);

        template <class type>
        static void recieve(type* buffer,int size,int source,int message_id) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
            MPI_Recv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        template <class type>
        static void send(type* buffer,int size,int source,int message_id) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
            MPI_Send(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD);
        }
        template <class type>
        static void allreduce(type* sendbuf,type* recvbuf,int size, Int op) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
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
        static void reduce(type* sendbuf,type* recvbuf,int size, Int op) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
            MPI_Op mpi_op;
            switch(op) {
                case OP_MAX: mpi_op = MPI_MAX; break;
                case OP_MIN: mpi_op = MPI_MIN; break;
                case OP_SUM: mpi_op = MPI_SUM; break;
                case OP_PROD: mpi_op = MPI_PROD; break;
            }
            MPI_Reduce(sendbuf,recvbuf,count,MPI_SCALAR,mpi_op,0,MPI_COMM_WORLD);
        }
        template <class type>
        static void irecieve(type* buffer,int size,int source,int message_id,void* request) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
            MPI_Irecv(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
        }
        template <class type>
        static void isend(type* buffer,int size,int source,int message_id,void* request) {
            int FLOAT_SIZE;
            MPI_Type_size(MPI_SCALAR, &FLOAT_SIZE);
            const int count = (size * sizeof(type) / FLOAT_SIZE);
            MPI_Isend(buffer,count,MPI_SCALAR,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
        }
        static void waitall(int count,void* request) {
            MPI_Waitall(count,(MPI_Request*)request,MPI_STATUS_IGNORE);
        }
};
#endif
