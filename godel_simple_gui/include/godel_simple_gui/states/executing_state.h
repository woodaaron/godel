#ifndef EXECUTING_STATE_H
#define EXECUTING_STATE_H

#include <actionlib/client/simple_action_client.h>
#include <godel_msgs/SelectMotionPlanAction.h>
#include "godel_simple_gui/gui_state.h"
#include <ros/ros.h>

namespace godel_simple_gui
{

class ExecutingState : public GuiState
{
  Q_OBJECT
public:
  // Constructor
  ExecutingState(const std::vector<std::string>& plans);

  // Entry and exit classes
  virtual void onStart(BlendingWidget& gui);
  virtual void onExit(BlendingWidget& gui);

  // Handlers for the fixed buttons
  virtual void onNext(BlendingWidget& gui);
  virtual void onBack(BlendingWidget& gui);
  virtual void onReset(BlendingWidget& gui);

protected:
  void executeAll();
  void executeOne(const std::string& plan);

private:
  bool execute_next_;
  std::vector<std::string> plan_names_;
  actionlib::SimpleActionClient<godel_msgs::SelectMotionPlanAction> select_motion_plan_action_client_;
  void selectMotionPlanDoneCallback(const actionlib::SimpleClientGoalState& state,
    const godel_msgs::SelectMotionPlanResultConstPtr& result);
};
}

#endif // EXECUTING_STATE_H
