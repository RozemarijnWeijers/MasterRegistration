#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"

class ClientIGT
{

  public:

  ClientIGT();
  ~ClientIGT();
  void ConnectToServer( char*, int );
  void DisconnectFromServer();
  int ReceiveImage();
  void SendImage();

  igtl::ClientSocket::Pointer socket;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  igtl::MessageHeader::Pointer headerMsg;
  igtl::TimeStamp::Pointer timeStamp;
  igtlUint32 sec;
  igtlUint32 nanosec;

};

