# include "usrClient.h"

Clients::Clients( char* host, int port )
{

  //Open a socket for th client
  socket = igtl::ClientSocket::New();

  //Connect to the server
  int r = socket->ConnectToServer( host, port );

  // Check the connection
  if ( (r) != 0 )
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    //exit(0);
  }

  this->headerMsg = igtl::MessageHeader::New();
  this->ts = igtl::TimeStamp::New();
  this->imgMsg = igtl::ImageMessage::New();

  std::cerr << "Client is connected to server: " << host << ":" << port << std::endl;

}

int Clients::ReceiveImage()
{

  // Initialize receive buffer
  this->headerMsg->InitPack();

  // Receive header message
  if (0 == this->socket->Receive(this->headerMsg->GetPackPointer(), this->headerMsg->GetPackSize()))
  {
    this->socket->CloseSocket();
    std::cerr<< " Exit 0 " <<std::endl;
    return 0;
  }

  // Unpack and deserialize the header
  headerMsg->Unpack();

  // Get time stamp
  headerMsg->GetTimeStamp(this->ts);
  //this->ts->GetTimeStamp(&sec, &nanosec);

  if (strcmp(this->headerMsg->GetDeviceType(), "IMAGE") == 0)
  {
    // Allocate memory for a message buffer to receive transform data
    this->imgMsg->SetMessageHeader(this->headerMsg);
    this->imgMsg->AllocatePack();
    // Receive transform data from the socket
    this->socket->Receive(imgMsg->GetPackBodyPointer(), this->imgMsg->GetPackBodySize());
    // Deserialize the transform data // If you want to skip CRC check, call Unpack() without argument.
    int c = this->imgMsg->Unpack(1);
    if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
      //image->SetParametersFromIGT();
      //image->IGTtoITKImage();
      std::cerr<<"Image received"<<std::endl;
      return 0;
    }
  }
  else
  {
    this->socket->Skip(this->headerMsg->GetBodySizeToRead(), 0);
  }

  return 1;

}

int Clients::SendImage()
{

  this->imgMsg->Pack();
  this->socket->Send( this->imgMsg->GetPackPointer(), this->imgMsg->GetPackSize() );

  return 1;

}
