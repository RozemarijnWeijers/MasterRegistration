# include "usrClientIGT.h"

ClientIGT::ClientIGT()
{

  //Open a socket for the client
  this->socket = igtl::ClientSocket::New();

  this->headerMsg = igtl::MessageHeader::New();
  this->timeStamp = igtl::TimeStamp::New();

}

ClientIGT::~ClientIGT()
{
}

void ClientIGT::ConnectToServer( char* serverHost, int serverPort )
{

  // Connect to Server and check the connection
  if ( this->socket->ConnectToServer( serverHost, serverPort ) != 0 )
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    return;
  }

  std::cerr << "ClientIGT is connected to server: " << serverHost << ":" << serverPort << std::endl;

  return;

}

void ClientIGT::DisconnectFromServer()
{

  this->socket->CloseSocket();

  return;

}

int ClientIGT::ReceiveImage()
{

  this->imgMsg = igtl::ImageMessage::New();

  // Initialize receive buffer
  this->headerMsg->InitPack();

  // Receive header message
  if (0 == this->socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize()))
  {
    this->socket->CloseSocket();
    return 1;
  }

  // Unpack and deserialize the header
  this->headerMsg->Unpack();

  // Get time stamp
  this->headerMsg->GetTimeStamp(timeStamp);
  this->timeStamp->GetTimeStamp(&sec, &nanosec);

  if (strcmp(this->headerMsg->GetDeviceType(), "IMAGE") == 0)
  {
    // Allocate memory for a message buffer to receive transform data
    this->imgMsg->SetMessageHeader(this->headerMsg);
    this->imgMsg->AllocatePack();
    // Receive transform data from the socket
    this->socket->Receive(this->imgMsg->GetPackBodyPointer(), this->imgMsg->GetPackBodySize());
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
    this->socket->Skip(headerMsg->GetBodySizeToRead(), 0);
  }

  return 1;

}

void ClientIGT::SendImage()
{

  this->imgMsg->Pack();
  this->socket->Send( imgMsg->GetPackPointer(), imgMsg->GetPackSize() );

  return;

}
