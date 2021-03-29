#include "smb_highlevel_controller/SmbHighlevelController.hpp"

namespace smb_highlevel_controller {

SmbHighlevelController::SmbHighlevelController(ros::NodeHandle& nodeHandle) :
  nodeHandle_(nodeHandle)
{
  if (!readParameters()) {
	ROS_ERROR("Could not read parameters.");
	ros::requestShutdown();
  }
  publisher_ = nodeHandle_.advertise<geometry_msgs::Twist>("/cmd_vel", 10);
  subscriber_ = nodeHandle_.subscribe(scan_topic, queue_size, &SmbHighlevelController::scanCallback, this);
  vis_pub = nodeHandle.advertise<visualization_msgs::Marker>( "visualization_marker", 0 );
  ROS_INFO("Successfully launched node.");
}

SmbHighlevelController::~SmbHighlevelController()
{
}

bool SmbHighlevelController::readParameters()
{
  if (!nodeHandle_.getParam("scan_topic", scan_topic)) return false;
  else if (!nodeHandle_.getParam("queue_size", queue_size)) return false;
  else if (!nodeHandle_.getParam("p_gain", p_gain)) return false;
  else if (!nodeHandle_.getParam("x_dot", x_dot)) return false;
  return true;
}

void SmbHighlevelController::scanCallback(const sensor_msgs::LaserScan& message)
{
   float smallestDistance = message.ranges[0];
   int smallest_i = 0;
   for (int i = 0; i < message.ranges.size(); i++){
	   if (smallestDistance > message.ranges[i]){
		   smallestDistance = message.ranges[i];
		   smallest_i = i;
	   }
   }
   ROS_INFO_STREAM("Smallest Distance: " << smallestDistance << std::endl);
   float angle_increment = message.angle_increment;
   ROS_INFO_STREAM("angle min: " << message.angle_min << std::endl);
   ROS_INFO_STREAM("angle increment: " << angle_increment << std::endl);
   float angle = message.angle_min + smallest_i * angle_increment;
   float x = smallestDistance * cos(angle);
   float y = smallestDistance * sin(angle);
   float z = 0;
   ROS_INFO_STREAM("x: " << x << std::endl << "y: " << y << std::endl);
   ROS_INFO_STREAM("angle: "<< angle << std::endl);
   float angle_dot = p_gain * (angle - 0);
   ROS_INFO_STREAM("angular velocity: "<< angle_dot << std::endl);

   float stop_distance = 2;
   if (x <= stop_distance){
   geometry_msgs::Twist stop_cmd;
   stop_cmd.linear.x = 0;
   stop_cmd.angular.z = 0;
   ROS_INFO_STREAM("cmd_vel: " << stop_cmd);
   publisher_.publish(stop_cmd);
   }
   else{
	   geometry_msgs::Twist msg;
	   msg.linear.x = x_dot;
	   msg.angular.z = angle_dot;
	   ROS_INFO_STREAM("cmd_vel: " << msg);
	   publisher_.publish(msg);
   }

   // visualization marker
   visualization_msgs::Marker marker;
   marker.header.frame_id = "base_link";
   marker.header.stamp = ros::Time();
   marker.ns = "my_namespace";
   marker.id = 0;
   marker.type = visualization_msgs::Marker::SPHERE;
   marker.action = visualization_msgs::Marker::ADD;
   marker.pose.position.x = x;
   marker.pose.position.y = y;
   marker.pose.position.z = z;
   marker.pose.orientation.x = 0.0;
   marker.pose.orientation.y = 0.0;
   marker.pose.orientation.z = 0.0;
   marker.pose.orientation.w = 1.0;
   marker.scale.x = 0.1;
   marker.scale.y = 0.1;
   marker.scale.z = 0.1;
   marker.color.a = 1.0; // Don't forget to set the alpha!
   marker.color.r = 0.0;
   marker.color.g = 1.0;
   marker.color.b = 0.0;
   //only if using a MESH_RESOURCE marker type:
   marker.mesh_resource = "package://pr2_description/meshes/base_v0/base.dae";
   vis_pub.publish(marker);
}

} /* namespace */
