import numpy as np
import quantities as pq
import neo


def add_channel_and_units_v2(block):
    sts = []
    for seg in block.segments:
        sts.extend(seg.spiketrains)
        try:
            del seg.events[0].annotations['signal']
        except:
            pass
    # collecting channel_indexes by channel_id
    channel_idxs = {}

    def get_channel_idx(channel_id):
        if channel_id not in channel_idxs:
            channel_idxs[channel_id] = neo.ChannelIndex(name='Channel {}'.format(channel_id),
    index=None,
    channel_id=channel_id)
            block.channel_indexes.append(channel_idxs[channel_id])
        return channel_idxs[channel_id]

    def get_unit(unit_id, channel_idx):
        unit = None
        for u in channel_idx.units:
            if u.annotations['unit_id'] == unit_id:
                unit = u

        if unit is None:
            unit = neo.Unit(unit_id=unit_id, channel_id=channel_idx.annotations['channel_id'])
            channel_idx.units.append(unit)

        return unit

    for st in sts:
        channel_id = st.annotations['channel_id']
        unit_id = st.annotations['unit_id']

        channel_idx = get_channel_idx(channel_id)
        unit = get_unit(unit_id, channel_idx)

        unit.spiketrains.append(st)