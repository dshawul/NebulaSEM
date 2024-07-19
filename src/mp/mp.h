#ifndef __MP_H
#define __MP_H

#include "mpi.h"
#include "types.h"

#if defined USE_DOUBLE
#   define SCALAR      double
#else
#   define SCALAR      float
#endif

#define  MY_PATH_MAX  512

// Define a traits class to map C++ types to MPI types
template <typename T>
struct mpi_traits;

template <>
struct mpi_traits<int> {
    static constexpr MPI_Datatype mpi_type = MPI_INT;
};
template <>
struct mpi_traits<unsigned> {
    static constexpr MPI_Datatype mpi_type = MPI_UNSIGNED;
};
template <>
struct mpi_traits<float> {
    static constexpr MPI_Datatype mpi_type = MPI_FLOAT;
};
template <>
struct mpi_traits<double> {
    static constexpr MPI_Datatype mpi_type = MPI_DOUBLE;
};

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

        template < typename type, typename typeb = SCALAR>
        static void recieve(type* buffer,int size,int source,int message_id) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Recv(buffer,count,mpi_traits<typeb>::mpi_type,source,message_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        template <typename type, typename typeb = SCALAR>
        static void send(type* buffer,int size,int source,int message_id) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Send(buffer,count,mpi_traits<typeb>::mpi_type,source,message_id,MPI_COMM_WORLD);
        }
        template <typename type, typename typeb = SCALAR>
        static void allreduce(type* sendbuf,type* recvbuf,int size, Int op) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Op mpi_op;
            switch(op) {
                case OP_MAX: mpi_op = MPI_MAX; break;
                case OP_MIN: mpi_op = MPI_MIN; break;
                case OP_SUM: mpi_op = MPI_SUM; break;
                case OP_PROD: mpi_op = MPI_PROD; break;
            }
            MPI_Allreduce(sendbuf,recvbuf,count,mpi_traits<typeb>::mpi_type,mpi_op,MPI_COMM_WORLD);
        }
        template <typename type, typename typeb = SCALAR>
        static void reduce(type* sendbuf,type* recvbuf,int size, Int op) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Op mpi_op;
            switch(op) {
                case OP_MAX: mpi_op = MPI_MAX; break;
                case OP_MIN: mpi_op = MPI_MIN; break;
                case OP_SUM: mpi_op = MPI_SUM; break;
                case OP_PROD: mpi_op = MPI_PROD; break;
            }
            MPI_Reduce(sendbuf,recvbuf,count,mpi_traits<typeb>::mpi_type,mpi_op,0,MPI_COMM_WORLD);
        }
        template <typename type, typename typeb = SCALAR>
        static void irecieve(type* buffer,int size,int source,int message_id,void* request) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Irecv(buffer,count,mpi_traits<typeb>::mpi_type,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
        }
        template <typename type, typename typeb = SCALAR>
        static void isend(type* buffer,int size,int source,int message_id,void* request) {
            const int count = (size * sizeof(type) / sizeof(typeb));
            MPI_Isend(buffer,count,mpi_traits<typeb>::mpi_type,source,message_id,MPI_COMM_WORLD,(MPI_Request*)request);
        }
        static void waitall(int count,void* request) {
            MPI_Waitall(count,(MPI_Request*)request,MPI_STATUS_IGNORE);
        }
};
#endif
