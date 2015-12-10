#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

class Clients
{

  public:

  Clients( char*, int );
  int ReceiveImage();
  int SendImage();

  igtl::ClientSocket::Pointer socket;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  igtl::MessageHeader::Pointer headerMsg;
  igtl::TimeStamp::Pointer ts;
  igtlUint32 sec;
  igtlUint32 nanosec;

};

