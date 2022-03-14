import numpy as np
import intervaltree

# Define bank of available fixed environments
environments_fixed = {
    "A":  [
        "SetReactionValue NOT 2.0",
        "SetReactionValue NAND 0.8",
        "SetReactionValue AND 2.0",
        "SetReactionValue ORN 0.8",
        "SetReactionValue OR 2.0",
        "SetReactionValue ANDN 0.8"
    ],
    "B": [
        "SetReactionValue NOT 0.8",
        "SetReactionValue NAND 2.0",
        "SetReactionValue AND 0.8",
        "SetReactionValue ORN 2.0",
        "SetReactionValue OR 0.8",
        "SetReactionValue ANDN 2.0"
    ]
}
fixed_env_order = ["A", "B"]
fixed_env_profiles = {"A":"101010", "B":"010101"}

tasks = ["NOT","NAND","AND","ORN","OR","ANDN"]
task_values = ["2.0", "0.8"]

def gen_const_env_events(total_updates, env):

    event_desc = []
    event_content = "\n"
    event_content += "\n".join([f"u begin {reaction}" for reaction in environments_fixed[env]])

    event_desc.append(
        {"profile": fixed_env_profiles[env], "label": env, "start_update": 0}
    )

    return {"event_content": event_content, "event_desc": event_desc}


def gen_cyclic_env_events(total_updates, cycle_period, seed, end_env_time=300):
    # Make this call repeatable for a given seed.
    rnd = np.random.default_rng(seed)

    event_content = ""
    event_desc = []

    total_updates += 1
    next_update = total_updates - end_env_time
    cur_state = 0
    # interval_tree = intervaltree.IntervalTree()
    # interval_tree[next_update:total_updates] = cur_state

    cur_env = fixed_env_order[cur_state]
    while (next_update > 0):
        # Write event content
        event_content += "\n"
        event_content += "\n".join([f"u {next_update} {reaction}" for reaction in environments_fixed[cur_env]])
        # Record change description
        event_desc.append(
            {"profile": fixed_env_profiles[cur_env], "label": cur_env, "start_update": next_update}
        )
        # prev_update = next_update
        next_update -= int(rnd.normal(loc=cycle_period, scale=cycle_period / 5))
        # interval_tree[next_update:prev_update] = cur_state

        cur_state = (cur_state + 1) % len(fixed_env_order)
        cur_env = fixed_env_order[cur_state]

    # Start with next state
    event_content += "\n"
    event_content += "\n".join([f"u begin {reaction}" for reaction in environments_fixed[cur_env]])
    # Record change description
    event_desc.append(
        {"profile": fixed_env_profiles[cur_env], "label": cur_env, "start_update": "0"}
    )

    return {"event_content": event_content, "event_desc": event_desc}

def gen_random_env_events(total_updates, chg_period, seed, end_env_time=300):
    rnd = np.random.default_rng(seed)

    event_content = ""
    event_desc = []

    total_updates += 1
    next_update = total_updates - end_env_time

    # Set final period to be env a
    cur_env = "A"
    event_content += "\n"
    event_content += "\n".join([f"u {next_update} {reaction}" for reaction in environments_fixed[cur_env]])
    event_desc.append(
        {"profile": fixed_env_profiles[cur_env], "label": fixed_env_profiles[cur_env], "start_update": next_update}
    )
    next_update -= int(rnd.normal(loc=chg_period, scale=chg_period / 5))

    # Now, switch to randomly changing
    while (next_update > 0):
        # Randomize the environment
        randomized_values = [float(rnd.choice(task_values)) for _ in tasks]
        env_profile = "".join([str(int(randomized_values[i] >= 1)) for i in range(len(tasks))])
        event_content += "\n"
        event_content += "\n".join([f"u {next_update} SetReactionValue {tasks[i]} {randomized_values[i]}" for i in range(len(tasks))])
        event_desc.append(
            {"profile": env_profile, "label": env_profile, "start_update": next_update}
        )

        next_update -= int(rnd.normal(loc=chg_period, scale=chg_period / 5))

    # Start with a random environment
    randomized_values = [float(rnd.choice(task_values)) for _ in tasks]
    env_profile = "".join([str(int(randomized_values[i] >= 1)) for i in range(len(tasks))])
    event_content += "\n"
    event_content += "\n".join([f"u begin SetReactionValue {tasks[i]} {randomized_values[i]}" for i in range(len(tasks))])
    event_desc.append(
        {"profile": env_profile, "label": env_profile, "start_update": "0"}
    )

    return {"event_content": event_content, "event_desc": event_desc}

if __name__ == "__main__":
    ret = gen_cyclic_env_events(10000, 30, 2)
    print(ret["event_content"])