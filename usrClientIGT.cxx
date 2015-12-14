# include "usrClientIGT.h"

ClientIGT::ClientIGT()
{

  //Open a socket for the client
  socket = igtl::ClientSocket::New();

  headerMsg = igtl::MessageHeader::New();
  timeStamp = igtl::TimeStamp::New();

}

ClientIGT::~ClientIGT()
{
}

void ClientIGT::ConnectToServer( char* serverHost, int serverPort )
{

  // Connect to Server and check the connection
  if ( socket->ConnectToServer( serverHost, serverPort ) != 0 )
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    return;
  }

  std::cerr << "ClientIGT is connected to server: " << serverHost << ":" << serverPort << std::endl;

  return;

}

void ClientIGT::DisconnectFromServer()
{

  socket->CloseSocket();

  return;

}

int ClientIGT::ReceiveImage()
{

  imgMsg = igtl::ImageMessage::New();

  // Initialize receive buffer
  headerMsg->InitPack();

  // Receive header message
  if (0 == socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize()))
  {
    socket->CloseSocket();
    return 1;
  }

  // Unpack and deserialize the header
  headerMsg->Unpack();

  // Get time stamp
  headerMsg->GetTimeStamp(timeStamp);
  timeStamp->GetTimeStamp(&sec, &nanosec);

  if (strcmp(headerMsg->GetDeviceType(), "IMAGE") == 0)
  {
    // Allocate memory for a message buffer to receive transform data
    imgMsg->SetMessageHeader(headerMsg);
    imgMsg->AllocatePack();
    // Receive transform data from the socket
    socket->Receive(imgMsg->GetPackBodyPointer(), imgMsg->GetPackBodySize());
    // Deserialize the transform data // If you want to skip CRC check, call Unpack() without argument.
    int c = imgMsg->Unpack(1);
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
    socket->Skip(headerMsg->GetBodySizeToRead(), 0);
  }

  return 1;

}

void ClientIGT::SendImage()
{

  imgMsg->Pack();
  socket->Send( imgMsg->GetPackPointer(), imgMsg->GetPackSize() );

  return;

}
