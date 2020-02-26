#!/usr/bin/env python
# -*- coding:utf-8-*-
# author: Liukewei time:2020/1/21 QQ:422209303 e-mail:Liukeweiaway@hotmail.com
# ----------------------------------------------------------------------------


def choiceThreshold(species, modification, threshold):
    threshold = str(threshold)
    if species == 'Human':
        if modification == 'm1A':
            if threshold == 'low':  # low
                M_threshold = 0.0237
            elif threshold == 'normal':  # normal
                M_threshold = 0.0368
            elif threshold == 'high':  # high
                M_threshold = 0.1071
        elif modification == 'm6A':
            if threshold == 'low':
                M_threshold = 0.2282
            elif threshold == 'normal':
                M_threshold = 0.2596
            elif threshold == 'high':
                M_threshold = 0.2882
        elif modification == 'm5C':
            if threshold == 'low':
                M_threshold = 0.7057
            elif threshold == 'normal':
                M_threshold = 0.8438
            elif threshold == 'high':
                M_threshold = 0.9052
        elif modification == 'pseudouridine':
            if threshold == 'low':
                M_threshold = 0.7314
            elif threshold == 'normal':
                M_threshold = 0.7683
            elif threshold == 'high':
                M_threshold = 0.8137
        elif modification == 'A-to-I':
            if threshold == 'low':
                M_threshold = 0.4484
            elif threshold == 'normal':
                M_threshold = 0.5286
            elif threshold == 'high':
                M_threshold = 0.6273

    elif species == 'Yeast':
        if modification == 'm1A':
            if threshold == 'low':
                M_threshold = 0.8039
            elif threshold == 'normal':
                M_threshold = 0.9419
            elif threshold == 'high':
                M_threshold = 0.9419
        elif modification == 'm6A':
            if threshold == 'low':
                M_threshold = 0.7314
            elif threshold == 'normal':
                M_threshold = 0.8130
            elif threshold == 'high':
                M_threshold = 0.8839
        elif modification == 'm5C':
            if threshold == 'low':
                M_threshold = 0.1429
            elif threshold == 'normal':
                M_threshold = 0.1838
            elif threshold == 'high':
                M_threshold = 0.7069
        elif modification == 'pseudouridine':
            if threshold == 'low':
                M_threshold = 0.6633
            elif threshold == 'normal':
                M_threshold = 0.7364
            elif threshold == 'high':
                M_threshold = 0.7706

    elif species == 'Mouse':
        if modification == 'm1A':
            if threshold == 'low':
                M_threshold = 0.0667
            elif threshold == 'normal':
                M_threshold = 0.0979
            elif threshold == 'high':
                M_threshold = 0.1225
        elif modification == 'm6A':
            if threshold == 'low':
                M_threshold = 0.5854
            elif threshold == 'normal':
                M_threshold = 0.6650
            elif threshold == 'high':
                M_threshold = 0.9196
        elif modification == 'm5C':
            if threshold == 'low':
                M_threshold = 0.2834
            elif threshold == 'normal':
                M_threshold = 0.9792
            elif threshold == 'high':
                M_threshold = 0.9792
        elif modification == 'pseudouridine':
            if threshold == 'low':
                M_threshold = 0.6723
            elif threshold == 'normal':
                M_threshold = 0.7058
            elif threshold == 'high':
                M_threshold = 0.7256
    return M_threshold
