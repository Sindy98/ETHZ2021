#pragma once

#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <geometry_msgs/Twist.h>
#include <visualization_msgs/Marker.h>
#include <stdlib.h>

namespace smb_highlevel_controller {

/*!
 * Class containing the Husky Highlevel Controller
 */
class SmbHighlevelController {
public:
	/*!
	 * Constructor.
	 */
	SmbHighlevelController(ros::NodeHandle& nodeHandle);

	/*!
	 * Destructor.
	 */
	virtual ~SmbHighlevelController();

private:
	ros::NodeHandle& nodeHandle_;
	ros::Subscriber subscriber_;
	ros::Publisher publisher_;
	ros::Publisher vis_pub;
	std::string scan_topic;
	std::string cmd_vel;
	int queue_size;
	float p_gain;
	float x_dot;
	bool readParameters();
	void scanCallback(const sensor_msgs::LaserScan& message);
};

} /* namespace */
