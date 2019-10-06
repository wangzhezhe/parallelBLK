//the rank 0 will gether the ip from all the server
//and then write it to current node

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <mpi.h>
#include <vector>

using namespace std;

//it is also ok to return map if we need support more complex using senarios
vector<string> split(const string &s, const string &seperator)
{
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size())
    {
        int flag = 0;
        while (i != s.size() && flag == 0)
        {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x])
                {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0)
        {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x])
                {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j)
        {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }
    return result;
}

string recordIPPort()
{
    string ipstr;
    int n;
    struct ifreq ifr;

    n = socket(AF_INET, SOCK_DGRAM, 0);
    //Type of address to retrieve - IPv4 IP address
    ifr.ifr_addr.sa_family = AF_INET;
    //Copy the interface name in the ifreq structure
    string INTERFACE = "eno1";
    strncpy(ifr.ifr_name, INTERFACE.data(), IFNAMSIZ - 1);
    ioctl(n, SIOCGIFADDR, &ifr);
    close(n);
    //display result
    char *ip = inet_ntoa(((struct sockaddr_in *)&ifr.ifr_addr)->sin_addr);

    ipstr = string(ip);
    return ipstr;
}

int main(int argc, char **argv)
{

    MPI_Init(NULL, NULL);

    // Get the rank of the process
    int rank;
    int procNum;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);

    string ip = recordIPPort();

    //maybe it is ok that only write the masternode's ip and provide the interface
    //of getting ipList for other clients

    if (rank == 0)
    {
        std::cout << "total process num: " << procNum << std::endl;
    }

    std::cout << "current rank is " << rank << " current ip is: " << ip << std::endl;

    //add extra position for the comma

    int sendLen = ip.length()+1;
    
    char*sendipStr=(char*)malloc(sizeof(char)*(sendLen));    
    sprintf(sendipStr, "%s,", ip.data()); 

    //std::cout << "check send ip: "<<string(sendipStr) << std::endl;

    int sendSize = sendLen* sizeof(char);

    int rcvLen = sendLen;
    int rcvSize = sendSize * procNum+10;
    char *rcvString = NULL;

    if (rank == 0)
    {
        std::cout << "sendSize: " << sendSize << ", rcvSize:" << rcvSize << std::endl;
        rcvString = (char *)malloc(rcvSize);
        {
            if (rcvString == NULL)
            {
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 0;
            }
        }
    }

    /*
    MPI_Gather(void* send_data,
    int send_count,
    MPI_Datatype send_datatype,
    void* recv_data,
    int recv_count,
    MPI_Datatype recv_datatype,
    int root,
    MPI_Comm communicator)
    */

    //attention, the recv part is the size of the buffer recieved from the each thread instead of all the str
    //refer to https://stackoverflow.com/questions/37993214/segmentation-fault-on-mpi-gather-with-2d-arrays
    int error_code = MPI_Gather(sendipStr, sendLen, MPI_CHAR, rcvString, rcvLen, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (error_code != MPI_SUCCESS)
    {
        std::cout << "error for rank " << rank << " get MPI_GatherError: " << error_code << std::endl;
    }
    //write to file for ip list json file if it is necessary
    //or expose the list by the rpc
    if (rank == 0)
    {
        //add termination for the last position
        rcvString[rcvLen*procNum-1]='\0';
        printf("rcv value: %s\n",rcvString);
        string list = string(rcvString);
        vector<string> ipList = split (list,",");
        
        std::cout << "check the ip list:" << std::endl;
        for(int i=0;i<ipList.size();i++){
            std::cout << ipList[i] << std::endl;
        }
        
        free(rcvString);
    }

    MPI_Finalize();
    return 0;
}